#!/bin/bash

# Test script for SHINE infant XGBoost models.

./run_child_job.R --debug \
    --parmfile=demo.rds \
    --jobname=test3 --jobid=interactive --array_index=1
