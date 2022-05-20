source( "xgboost.R" );
library( glue );

select <- dplyr::select;

spread_relab <- function( data, exclude=c( "Clade","NCBI_ID","Cov","n_Reads","rank" ), RelAb.fill=0, tfun=identity )
    pivot_wider( data, names_from=name, values_from=RelAb, values_fn=tfun, values_fill=list( RelAb=RelAb.fill ),
               id_cols=-one_of( exclude ) );

# spread_relab <- function( data, id_cols=-one_of( "Clade","NCBI_ID","Cov","n_Reads","rank" ) )
cummulative_n_distinct <- function( x )
{
    map( 1:length( x ), partial( head, x ) ) %>% map_int( n_distinct );
}

# This pulls levels of enum columns from H2O frame.
filter_varimp <- function( frameId, models, features, model.score.quantile=0.05, importance.threshold=0.95, n_vars_min=20 )
{
    hf1 <- h2o.getFrame( frameId );
    enum1 <- h2o.getTypes( hf1 ) %>% unlist() %>% { which( . == "enum" ); }

#   Construct a map from {variable}.{level} to {variable}.
    variable_map <-
        as_tibble( hf1[ , enum1 ] ) %>%
        mutate( across( everything(), fct_explicit_na, na_level = "missing(NA)" ) ) %>%
        map( levels ) %>% enframe( value = "level" ) %>% unnest( level ) %>%
        mutate( variable = glue( "{name}.{level}" ) );

    # Models are BayesOptimization iterates;  They all have the same set of features.
#   features <- pluck( models, "x", 1 );

    feat1 <-
        enframe( features, name = "feature_group", value = "variable" ) %>%
        mutate( across( feature_group, compose( fct_inorder, factor ) ) );

#   Select the best models.  Lower scores indicate better models.
    filter( models, Value <= quantile( Value, probs=model.score.quantile ) ) %>%
    unnest( varimp ) %>%
    left_join( variable_map ) %>%
#   Strip factor level from variable column.
    mutate( variable = ifelse( is.na( name ), variable, name ) ) %>%
    left_join( feat1 ) %>%
    group_by( model_id, feature_group ) %>% arrange( desc( scaled_importance ) ) %>%
#   Apply importance threshold filter by feature group.
    mutate( p_group = cumsum( scaled_importance ) / sum( scaled_importance ) ) %>%
    select( variable, scaled_importance, p_group ) %>%
    ungroup( feature_group ) %>% arrange( desc( scaled_importance ) ) %>%
    mutate(
#       Apply importance threshold filter overall.
        p = cumsum( scaled_importance ) / sum( scaled_importance ),
        k = cummulative_n_distinct( variable )
    ) %>%
#   Retain each feature above importance threshold by either feature group or overall.
    filter( ( p < importance.threshold ) | ( p_group < importance.threshold ) | ( k <= n_vars_min ) ) %>%
    ungroup() %>%
    distinct( feature_group, variable ) %>%
    arrange( feature_group ) %>% deframe();
}

get_h2o_metrics <- function( x )
{
    c( mse=h2o.mse, mae=h2o.mae, r2=h2o.r2, rmse=h2o.rmse, rmsle=h2o.rmsle ) %>%
    map( possibly, otherwise=NA ) %>% map( invoke, list( x ), xval=TRUE ) %>% as.list() %>% as_tibble();
}

three_stage_bayes_grid <- function( frameId, response="laz", features=NULL, init_points=10, n_iter=50, seed=42,
    final_run_time=1200, nfolds=c( 10,0,0 ), h2o_metric = NA, backend=c( "auto","gpu","cpu" ), .path=NULL )
{
    backend <- match.arg( backend );
    call.parms <- as.list( environment() );

    # This returns total CPU time in milliseconds for all cross-validation models.
    get_cv_runtime <- function( model_id, nfolds )
    {
        tibble(
            model_id = paste( model_id, "cv", 1:nfolds, sep='_' ),
            model = map( model_id, h2o.getModel ) %>% map( pluck, "model" ),
            run_time = map_dbl( model, pluck, "run_time" )
        ) %>% with( sum( run_time ) );
    }

    get_fold_count <- function( model )
    {
        # For classification models, which use fold_column '.folds', nfolds is zero.  
        # Accessing fold_column is also problematic.  Instead, count the number of cross-validation models.
        pluck( model, "model", "cross_validation_models" ) %>% length();
    }

    get_holdout_predictions <- function( model )
    {
        folds <- h2o.cross_validation_fold_assignment( model );
        # sum will set the column name to 'sum'.
        predict <- h2o.cross_validation_predictions( model ) %>% invoke( h2o.cbind, . ) %>%
            h2o.sum( axis=1, return_frame = TRUE );
        names( predict ) <- "predict";
        h2o.cbind( predict, folds ) %>% as_tibble();
    }

    summarize_models <- function( dt )
    {
#       Nest x, the vector of model features, before converting to tibble.
#       get_parameters <- function( model ) within( model@allparameters, { x <- list( x ) } ) %>% as_tibble();
#       There may be more than one list/vector column, e.g. x, fold_column.
        get_parameters <- function( model ) pluck( model, "allparameters" ) %>% map_if( negate( is.scalar ), list ) %>% as_tibble();

        as_tibble( dt ) %>%
        mutate(
            model = map( model_id, h2o.getModel ),
            allparameters = map( model, get_parameters ),
#           dt contains some of the parameters;  avoid duplicate names.
            across( allparameters, map, select, -any_of( names( . ) ) ),
            varimp = map( model, h2o.varimp ),
            holdout_predictions = map( model, get_holdout_predictions )
        ) %>% unnest( allparameters ) %>%
        mutate(
            nfolds = map_int( model, get_fold_count ),
            run_time = map2_dbl( model_id, nfolds, get_cv_runtime )
        );
    }
                 
    save_models <- function( ids, path )
    {
        res <- tibble(
            modelId = ids,
            model = map( modelId, h2o.getModel ),
#           model = llply( modelId, h2o.getModel, .parallel=F ),
            allparameters = map( model, slot, "allparameters" ),
            frame = map( allparameters, getElement, "training_frame" ) %>% map_chr( as.character ),
            score = map_dbl( model, h2o.mean_residual_deviance, xval=T ),
            varimp = map( model, h2o.varimp ),
            ntrees = map_int( allparameters, pluck, "ntrees" ),
            max_depth = map_int( allparameters, pluck, "max_depth" ),
            nfolds = map_int( model, get_fold_count ),
#           nfolds = map_int( allparameters, pluck, "nfolds" ),
            run_time = map2_dbl( modelId, nfolds, get_cv_runtime ),
            holdout_prediction = get_holdout_predictions( model )
        );

        select( res, -model ) %>% return();

#       Don't write mojo.
        walk( res$model, h2o.save_mojo, path );

        return( res );
    }

    if ( is.null( .path ) )
        .path <- "./results";
    
    rand_str <- sample( LETTERS, 5, replace=T ) %>% paste( collapse='' );
    out.path <- glue( "{.path}/child_{response}_xgboost_{t}_{rand_str}", t=format( Sys.time(), "%F_%Hh%M" ) );
    dir.create( out.path );

    # All elements of features list have the same set of id cols.  Remove the names.
    id_cols <- split( features, names( features ) ) %>% pluck( "id" ) %>% set_names( NULL );

    # Hyperparameter search with Bayesian Optimization: Stage 1
    
    r1 <- build_xgboost_model( frameId, response=response, type="grid", nfolds=nfolds[1], seed=seed,
            init_points=init_points, n_iter=n_iter, backend=backend, features=features, h2o_metric = h2o_metric );

    m1 <- summarize_models( r1$model$bayes$result$History );

    # Stage 2
    features2 <- filter_varimp( frameId, m1, features );

    # Hyperparameter search with Bayesian Optimization: Stage 2
    r2 <- build_xgboost_model( frameId, response=response, type="grid", nfolds=nfolds[2], seed=seed,
                init_points=init_points, n_iter=n_iter, features=features2, h2o_metric = h2o_metric );

    m2 <- summarize_models( r2$model$bayes$result$History );

    # Hyperparameter search with Bayesian Optimization: Stage 3
    features3 <- filter_varimp( frameId, m2, features );

    # Estimate ntrees required for final stage of desired duration.
    fit1 <- lm( run_time ~ ntrees + max_depth + I( ntrees*max_depth ), data=m2 );
    fit1.coef <- set_names( fit1$coefficients, letters[1:4] ) %>% as.list();
   
    best_parms <- pluck( r2,"model","bayes","result","best_parms" ) %>%
        with(
            {
                c1 <- ntrees * learn_rate;
                # multiply by 1000 to convert to milliseconds.  Assume average load of 20% of available threads.
                # Limit number of trees because actual runtimes are too long.
                n <- with( fit1.coef,
                           ( 1000 * final_run_time * ( 0.2 * nthreads ) - a - c*max_depth ) / ( b + d*max_depth )
                     ) %>% min( 10000 );
                within( .,
                    {
                        # Increase number of trees using linear model with a minimum bound.
                        ntrees <- as.integer( n );

                        # Increase ntrees, holding, ntrees * learn_rate, constant.
                        learn_rate <- c1/ntrees;
                        names( ntrees ) <- NULL;
                        names( learn_rate ) <- NULL;
                    }
                );
            }
        );
    best_parms$log10_ntrees <- NULL;
    best_parms$log10_learn_rate <- NULL;

    r3 <- rlang::dots_list( frameId=frameId, response=response, type="single", !!! best_parms, nfolds=nfolds[3], seed=seed,
                features=features3, .homonyms = "last" ) %>%
          invoke( build_xgboost_model, . );

    m3 <- pluck( r3,"model","res","result","model_id" ) %>% tibble( model_id = . ) %>% summarize_models();

    res <- tibble(
        stage = 1:3,
        r = list( r1, r2, r3 ),
        m = list( m1, m2, m3 ),
        features = list( features, features2, features3 )
    );
    return( res );
}

