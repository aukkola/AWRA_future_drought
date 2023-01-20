#!/bin/bash

#PBS -P w97
#PBS -q normal
#PBS -l storage=gdata/er4+scratch/w97+gdata/hh5+gdata/wj02+gdata/eg3+gdata/w97
#PBS -N drought_metric_job_BIAS-METHOD_MODEL_SCENARIO_VARIABLE_REGION 
#PBS -l walltime=03:00:00
#PBS -l ncpus=1
#PBS -l mem=70gb
#PBS -j oe


module use /g/data/hh5/public/modules
module load conda/analysis3

python3 /g/data/w97/amu561/Steven_CABLE_runs/scripts/drought_scripts/template_drought_metric_complete_24_month.py "BIAS-METHOD" "MODEL" "SCENARIO" "VARIABLE"
