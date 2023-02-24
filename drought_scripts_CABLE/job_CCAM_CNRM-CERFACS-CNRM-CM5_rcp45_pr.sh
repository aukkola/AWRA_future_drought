#!/bin/bash

#PBS -P w97
#PBS -q express
#PBS -l storage=gdata/er4+scratch/w97+gdata/hh5+gdata/wj02+gdata/eg3+gdata/w97
#PBS -N drought_metric_job_CCAM_CNRM-CERFACS-CNRM-CM5_rcp45_pr_REGION 
#PBS -l walltime=06:00:00
#PBS -l ncpus=24
#PBS -l mem=10gb
#PBS -j oe


module use /g/data/hh5/public/modules
module load conda/analysis3

mkdir -p test

python3 /g/data/w97/amu561/Steven_CABLE_runs/drought_scripts/template_drought_metric_complete_3_month.py "CCAM" "CNRM-CERFACS-CNRM-CM5" "rcp45" "pr"
