#!/usr/bin/env Rscript

.libPaths();

# usage: run_child_job [--debug] --parmfile=<parameter-file> --array_index=<job_array_index> --jobid=<jobid> --feature_set=<feature_set> --results_dir=<results_dir> --temp_dir=<temp_dir>
#       --response=<response>

suppressPackageStartupMessages(
    {
        library( R.utils );
        library( magrittr );
        library( tidyverse );
        library( assertthat );
        library( glue );
        library( h2o );
    }
);

source( "three_stage_booster.R" );

safely_as_integer <- function( x, default=NA )
    with( safely( as.integer )( x ), if ( is.null( result ) || is.na( result ) || ( length( result ) == 0 ) ) default else result );

concat_parms <- function( parms )
{
    enframe( parms ) %>% pmap( function( name, value ) glue( "{name}={value}" ) ) %>% paste( collapse=':' );
}

if ( ! exists( "h2o_server" ) )
    h2o_server <- list( ip = "localhost", port = -1 );

#start_h2o_server <- function( temp_dir )
start_h2o_server <- function()
{
        # Is H2O server instance active?
        # When run as an array job, result of tempdir() will contain special characters '[]'.
        # Using file.symlink obviates the need for backslash-enquoting.
        # Home is mounted read-only on job nodes, so create symlink in scratch.
#       h2o_tmp <- tempfile( pattern="ice_root.", tmpdir=glue( "{temp_dir}/h2o" ) );
#       file.symlink( tempdir(), h2o_tmp );
        find_available_port <- "import socket; s=socket.socket(); s.bind(('',0)); print( s.getsockname()[1] ); s.close()";
        cmd <- glue( "python -c \"{find_available_port}\"" );
        print( cmd );
        port <- system( cmd, intern=TRUE ) %>% as.integer();
        # Set bind_to_localhost=FALSE so H2O server port binds to all network interfaces.
#       h2o.init( port=port, ice_root=h2o_tmp, bind_to_localhost=FALSE );
        h2o.init( port=port, bind_to_localhost=FALSE );
        h2o_server$port <<- port;    
}

.Last <- function()
{
    # Shut down H2O server instance.
    if ( h2o_server$port >= 0 ) {
        with( h2o_server, h2o.connect( ip, port ) );
        h2o.shutdown();
    }
}

main <- function( argv )
{
    assert_that( has_name( argv, "parmfile" ), msg="Please specify parameter file with '--parmfile=' option." );
    assert_that( has_name( argv, "array_index" ), msg="Please specify array job index with '--array_index=' option." );
    assert_that( has_name( argv, "jobid" ), msg="Please specify job id with '--jobid=' option." );
    
    if ( is.null( argv$results_dir ) )
        argv$results_dir <- ".";

    argv$array_index %<>% safely_as_integer();

    # Connect to existing H2O server or start a new one.
    if ( h2o_server$port < 0 ) start_h2o_server() else with( h2o_server, h2o.connect( ip, port ) );

    p1 <- read_rds( argv$parmfile ) %>% pluck( "data" ) %>% map( argv$array_index );

    nest_nonscalars <- function( x )
        ifelse( is.null( x ) || ( length( x ) > 1 ), list( x ), x )

    hf1 <- spread_relab( p1$data ) %>% as.h2o();
    hf1_id <- h2o.getId( hf1 );
    
    print( glue( "debug={argv$debug}" ) );

    p0 <-
        if ( argv$debug )
            list( nfolds = c( 5,5,5 ), n_iter = 2, final_run_time = 100 )
        else list( nfolds = c( 5,0,0 ) ); 

    model_helper <- function( frameId, response, features, h2o_metric = "h2o.mse", ... )
    {
        r1 <- invoke( safely( three_stage_bayes_grid ),
                        rlang::dots_list( frameId, response, features, h2o_metric = h2o_metric, !!! p0 ) );
        map( r1, list ) %>% as_tibble();
    }

    r1 <- p1$parms %>%
        mutate(
            argv = list( argv ),
            frameId = hf1_id
        ) %>%
        mutate(
            models = pmap( ., model_helper )
        ) %>%
        unnest( models );

    # exitval is TRUE iff an error has occurred.
    
    exitval <- map_lgl( r1$error, negate( is.null ) ) %>% any();
    if ( exitval )
        print( r1$error )
    else {
        # Add random string to avoid collision when multiple jobs start at the same time.
        rand_str <- sample( LETTERS, 5, replace=T ) %>% paste( collapse='' );
        now <- format( Sys.time(), "%F_%Hh%M" );
        # Remove parms that won't be included in file name.
#       p2 <- parms[ ! names( parms ) %in% c( "data","features","parmfile" ) ];
        dir <- glue( "{argv$results_dir}/{argv$jobid}_{now}_{rand_str}" );
        print( glue( "Results will be written to directory '{dir}'." ) );

        # The dir will be created by one array job.  Don't try to create it multiple times.
        if ( !dir.exists( dir ) )
            dir.create( dir ); 

#       Save only the final stage model.
        save_all_mojo <- function( x, .path )
            tail( x,1 ) %>% unnest( m ) %>% pull( model ) %>% map_chr( h2o.save_mojo, path = .path );

        # Save each model in MOJO format.
        # Rename to avoid duplicate column names while unnesting.
        r1 %<>% mutate( mojo = map( result, save_all_mojo, dir ) );

        # Save nested tibble of model results.
        fname <- glue( "{dir}/child_model_results.rds.gz" );
        print( glue( "Results will be written to file '{fname}'." ) );
#       fname <- with( p2, glue( "{dir}/child_", concat_parms( p2 ), ".rds.gz" ) );
        write_rds( r1, fname, compress = "gz" );
        print( glue( "Model results saved to '{fname}'." ) );

        # Save H2O frame.
        h2o.exportFile( hf1, glue( "{dir}/{hf1_id}.csv.gz" ), compression = "gzip" );
    }

   h2o.shutdown( prompt=FALSE );

    return( exitval );

}

if ( !interactive() ) {
    .defaults <- list( debug = FALSE );
    quit( status=main( commandArgs( trailingOnly = TRUE, asValues = TRUE, defaults = .defaults ) ) );
}
