#!/bin/bash

#PBS -P oq98
#PBS -q normal
#PBS -j oe
#PBS -l storage=gdata/er4+scratch/w97+gdata/hh5+gdata/wj02+gdata/w97+gdata/oq98
#PBS -N drought_metric_job_BIAS-METHOD_MODEL_SCENARIO_VARIABLE_REGION 
#PBS -l walltime=00:10:00
#PBS -l ncpus=1
#PBS -l mem=1gb

PBS_PATH=/g/data/oq98/amu561/Steven_CABLE_runs/scripts/drought_scripts/pbs_jobs

mkdir -p ${PBS_PATH}

cd /g/data/oq98/amu561/Steven_CABLE_runs/scripts/drought_scripts

#for bcm in MRNBC; do
for bcm in QME MRNBC CCAM ISIMIP2b; do
    #for model in CNRM-CERFACS-CNRM-CM5; do
    for model in CSIRO-BOM-ACCESS1-0 CNRM-CERFACS-CNRM-CM5 NOAA-GFDL-GFDL-ESM2M MIROC-MIROC5; do
        #for scenario in rcp85; do
        for scenario in rcp45 rcp85; do
            #for variable in qtot s0; do
            for variable in qtot s0 pr; do
                              
		            outdir="/g/data/w97/amu561/Steven_CABLE_runs/drought_metrics/3-month/${bcm}/${model}/"  
                mkdir -p $outdir
		            out_file=/g/data/w97/amu561/Steven_CABLE_runs/drought_metrics/3-month/${bcm}/${model}/drought_metrics_${bcm}_${model}_${variable}_${scenario}_3.nc 
                     
                if [ ! -f ${out_file} ]; then

                echo job_${bcm}_${model}_${scenario}_${variable}
                cp job_template_3.sh ${PBS_PATH}/job_${bcm}_${model}_${scenario}_${variable}.sh
                wait
                sed -i "s/BIAS-METHOD/${bcm}/g" ${PBS_PATH}/job_${bcm}_${model}_${scenario}_${variable}.sh 
                wait 
                sed -i "s/MODEL/${model}/g" ${PBS_PATH}/job_${bcm}_${model}_${scenario}_${variable}.sh 
                wait
                sed -i "s/SCENARIO/${scenario}/g" ${PBS_PATH}/job_${bcm}_${model}_${scenario}_${variable}.sh 
                wait
                sed -i "s/VARIABLE/${variable}/g" ${PBS_PATH}/job_${bcm}_${model}_${scenario}_${variable}.sh 
                wait
                qsub ${PBS_PATH}/job_${bcm}_${model}_${scenario}_${variable}.sh
                    
                fi    
            done
        done
    done
done


