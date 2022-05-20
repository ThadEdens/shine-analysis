library( magrittr );
library( tidyverse );
library( assertthat );

library( rBayesianOptimization );
library( h2o );

nthreads <- 36;

options( "h2o.use.data.table"=TRUE );

# Argument features should be a list of character vectors of column names.  One model will be run for each list element.
build_xgboost_model <- function( frameId, response=c( "laz","waz","age" ), RelAb.fill=0, type=c( "single","grid" ),
    features=NULL, ... )
{
    .type <- match.arg( type );

#   Found at, https://stackoverflow.com/questions/11885207/get-all-parameters-as-list
#   argv <-  c( as.list( environment() ), list( ... ) );
    argv <- within( as.list( environment() ), dots <- list( ... ) );

    f1 <- switch( .type, single=xgb_model, grid=tune_grid_bayes )

    model_helper <- function( method )
    {
        function( frameId, response, features, ... )
        {
#           Don't pass on ... but do pass on ... of grandparent function.
            argv <-  c( as.list( environment() ), argv$dots );
            invoke( method, argv );
        }
    }

    mod1 <- model_helper( f1 )( frameId, response, features );

    res <- list(
                frameId = frameId,
                response = response,
                features = features,
                type = type,
                model = mod1
           );

    return( res );
}

xgb_model <- function( frameId, response, features, nfolds=0, seed=42, fold_column=NULL, ... )
{
#   Append dots to formal arguments, keeping last instance of each variable.
    argv <- rlang::dots_list( !!! as.list( environment() ), !!! rlang::list2( ... ), .homonyms="last" );

    hf1 <- h2o.getFrame( frameId );
    nobs <- h2o.nrow( hf1 );

    argv %<>% within(
        {
            if ( is.null( fold_column ) && ".folds" %in% h2o.names( hf1 ) )
                fold_column <- ".folds";

            if ( is.null( fold_column ) )
            {
                fold_assignment <- if ( nfolds==0 ) "Modulo" else "AUTO";
                if ( nfolds==0 )
                    nfolds <- nobs;
            } else {
                nfolds <- NULL;
                fold_assignment <- NULL;
            }

            training_frame <- frameId;
            y <- response;
            x <- setdiff( features, response );
            keep_cross_validation_predictions <- TRUE;
            keep_cross_validation_fold_assignment <- TRUE;

            frameId <- NULL;
            features <- NULL;
            response <- NULL;
        }
    ) %>% Filter( Negate( is.null ), . );
                     
    res <- invoke( safely( h2o.xgboost ), argv );
    list( argv=argv, res=res ) %>% return();
}

# The function, BayesianOptimization, exits without error ( just reporting stopping time ) when init_points is small.
# This does not appear to be documented.  I suspect that init_points must not be much smaller then length( grid_bounds ).
tune_grid_bayes <- function( frameId, response, features,
    nfolds = 0, init_grid_dt = NULL, init_points = 10, n_iter = 50, seed = 42, n_cv_repeats = 1,
    h2o_metric = c( NA, "h2o.mse", "h2o.logloss", "h2o.auc", "h2o.aucpr" ),
    fold_column = NULL, ... )
{
#   Add optional arguments to current environment.  This will update existing variables in current environment
#   with values from dots.
    list2env( list( ... ), environment() );

    hf1 <- h2o.getFrame( frameId );
    nobs <- h2o.nrow( hf1 );

#   If frame has a fold column, use it in the model.  This is passed to balance categorical response across folds.
    if ( is.null( fold_column ) && ".folds" %in% h2o.names( hf1 ) )
        fold_column <- ".folds";
    if ( is.null( fold_column ) )
    {
        fold_assignment <- if ( nfolds==0 ) "Modulo" else "AUTO";
        if ( nfolds==0 ) { nfolds <- nobs };
    } else {
        nfolds <- NULL;
        fold_assignment <- NULL;
    }

    response.type <- h2o.getTypes( hf1[ response ] )[[1]];

    model_ids <- c();

#   AUC, resp. MSE, is larger, resp. smaller, for better models.
#   Multinomial xgboost model does not have AUC whereas both binomial and multinomial have logloss.
    n_categories <- as.vector( hf1[ response ] ) %>% factor() %>% levels() %>% length();

    change_sign <- function( f ) compose( partial( `*`, -1.0 ), f );
    on_xval <- function( f ) partial( f, xval = TRUE );

    h2o_metric <- match.arg( h2o_metric );

    if ( is.na( h2o_metric ) )
        h2o_metric <- if ( response.type == "enum" ) "h2o.auc" else "h2o.mse";

    scoring_fn <- partial( get( h2o_metric ), xval = TRUE );

    if ( h2o_metric %in% c( "h2o.mse","h2o.logloss" ) )
        scoring_fn <- change_sign( scoring_fn );
    
    loss_bayes <- function( max_depth, log10_ntrees, log10_learn_rate, sample_rate, col_sample_rate, col_sample_rate_per_tree,
                min_child_weight )
    {
        xgb_args = list(
            x = features[ names( features ) != "id" ] %>% setdiff( response ),
            y = response,
            training_frame = frameId,
#           max_runtime_secs = 1,   # for testing
            fold_column = fold_column,
            nfolds = nfolds,
            fold_assignment = fold_assignment,
            keep_cross_validation_predictions = TRUE,
            keep_cross_validation_fold_assignment = TRUE,
            seed = seed,
            max_depth = max_depth, # rlang::int(max_depth),
            ntrees = round( 10.0 ** log10_ntrees, 0 ) %>% as.integer(),
#           min_rows = int(min_rows),
            learn_rate = 10.0 ** log10_learn_rate,
            sample_rate = sample_rate,
            col_sample_rate = col_sample_rate,
            col_sample_rate_per_tree = col_sample_rate_per_tree,
            min_child_weight = min_child_weight
        );
        fit <- invoke( h2o.xgboost, xgb_args );
        
#       Save model id.
        model_ids <<- c( model_ids, fit@model_id );

        list(
            Score = scoring_fn( fit ),
            Pred = h2o.cross_validation_predictions( fit )
        ) %>% return();
    }

    # When h2o.xgboost is called with max_depth=1, each tree will be a single split.  Confidence intervals on ALE plots
    # for such a model would be very small.
    grid_bounds <- list(
        max_depth = c( 2L,5L ), 
        log10_ntrees = c( 1.5, 2.5 ),
        log10_learn_rate = c( log10( 0.1 ), log10( 0.3 ) ),
        sample_rate = c( 0.8,1 ),
        col_sample_rate = c( 0.3,1 ),
        col_sample_rate_per_tree = c( 0.5, 1.0 ),       # multiplicative with col_sample_rate
        min_child_weight = c( 1L,3L )
    );

    k1 <- sum( init_points, nrow( init_grid_dt ) );
    assert_that( k1 >= length( grid_bounds ),
        msg = glue( "tune_grid_bayes: Stopping because init_points + nrow( init_grid_dt ), {k1},",
            " is less than length( grid_bounds ), { length( grid_bounds ) }." )
    );

    res <- safely( BayesianOptimization )( loss_bayes, bounds=grid_bounds, init_grid_dt=init_grid_dt,
                init_points=init_points, n_iter=n_iter );

    if ( !is.null( res$result ) )
    {
        res$result$best_parms <-
            within( as.list( res$result$Best_Par ),
                {
                    ntrees <- as.integer( round( 10 ** log10_ntrees ) );
                    learn_rate <- 10 ** log10_learn_rate;
                }
            );
        res$result$History$model_id <- model_ids;
        res$result$h2o_metric <- h2o_metric;
    }

    list( bounds=grid_bounds, bayes=res ) %>% return();
}
