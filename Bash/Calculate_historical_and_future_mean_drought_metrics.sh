




#Set path
path="/g/data/w97/amu561/Steven_CABLE_runs/"

#Find scales
scale=`ls $path/drought_metrics`


#Set experiments
declare -a exp=("historical" "rcp45" "rcp85")

#Set variables
declare -a vars=("pr", "qtot", "sm")

#Set metrics
declare -a metrics=("duration", "rel_intensity")



#Loop through scales
for s in $scale
do
  
  #Find bias correction methods
  bc_methods=`ls $path/drought_metrics/$s`
  
  #Loop through bc methods
  for b in $bc_methods
  do
    
    #Find GCMs
    gcms=`ls $path/drought_metrics/$s/$b`

    #Loop through variables
    for v in ${vars[@]}
    do
      
      #Loop through experiments
      for e in ${exp[@]}
      do
        
        #Find data file
        infile=`ls $path/drought_metrics/$s/$b/`
        
        #Set time period
        if [ $e = "historical" ]; then
          start_yr=1970
          end_yr=2005
        else
          start_yr=2064
          end_yr=2099
        fi
          
        #Loop through metrics
        for m in ${metrics[@]}
        do
          
          outfile= 
          
          #Select metric and time period and calculate mean
          cdo timmean -selyear,$start_yr/$endyr -selvar,$m $infile $outfile 
        
        
        done #metrics
        
      done #experiments
    
    done #vars
    
  done #bc_methods

done #scale
