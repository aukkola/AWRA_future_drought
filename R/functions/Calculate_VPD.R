#Calculate Qair from AWRA method based on mean and min daily temperature

#Equations 35-36 in: 
#The Australian Landscape Water Balance model (AWRA-L v6)
#Technical Description of the Australian Water Resources Assessment Landscape model version 6

#These give VPD, convert to Qair using FluxnetLSM method

#Assume constant air pressure of 100,000 Pa (see Mengyuan's GW paper)

#AWRA functions are modified to use tmin for both saturated and actual
#vapour pressure to obtain "9am" VPD value for weather generator
#and tmax to obtain "3pm" value

calculate_vpd <- function(tmax, tmin, air_pressure=100000) {
  
  
  #AWRA inputs are in Kelvin, convert to C here
  #tair <- tair - 273.15
  tmax  <- tmax - 273.15
  tmin  <- tmin - 273.15
  
  #Saturation vapour pressure (eq. 35), in Pa
  sat_pressure <- 610.8 * exp((17.27 * tmax) / (237.3 + tmax)) #610.8 * exp((17.27 * tmean) / (237.3 + tmean)) 
  
  
  #Actual vapour pressure (eq. 36), in Pa
  act_pressure <- 610.8 * exp((17.27 * tmin) / (237.3 + tmin)) #610.8 * exp((17.27 * tmin) / (237.3 + tmin))
  
  
  #Vapour pressure deficit (Pa)
  vpd <- sat_pressure - act_pressure
  
  #Convert from Pa to Hpa
  vpd_hPa <- vpd * 0.01
    
  
  # #Relative humidity (%)
  # relHum <- 100 * (1 - (vpd / sat_pressure))
  # 
  # 
  # # Then specific humidity at saturation:
  # ws <- 0.622 * sat_pressure / (air_pressure - sat_pressure)
  # 
  # # Then specific humidity:
  # specHum <- (relHum/100) * ws
  # 
  # 
  
  return(vpd_hPa)
  
}




