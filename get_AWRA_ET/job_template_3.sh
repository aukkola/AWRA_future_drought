#!/bin/bash

#PBS -P oq98
#PBS -q normal
#PBS -l storage=gdata/er4+scratch/w97+gdata/hh5+gdata/wj02+gdata/eg3+gdata/w97+gdata/oq98
#PBS -N drought_metric_job_BIAS-METHOD_MODEL_SCENARIO_VARIABLE_REGION 
#PBS -l walltime=02:00:00
#PBS -l ncpus=1
#PBS -l mem=60gb
#PBS -j oe


module use /g/data/hh5/public/modules
module load conda/analysis3

python3 /g/data/oq98/amu561/Steven_CABLE_runs/scripts/get_AWRA_ET/AWRA_future_ET.py "BIAS-METHOD" "MODEL" "SCENARIO" "VARIABLE"
