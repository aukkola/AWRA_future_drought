

#CHECK cos and sin functions if not working !!!!!

#SWdown in MJ/m2/day (Kd in Section 3.2.3, first equation)

#Vapour pressure in Pa (ea in Section 3.2.3, first equation)

#Temperature input in Kelvin (equations modified accordingly by 
#removing conversion to K (+ 273.15))

calculate_lwdown <- function(tmean_K, tmin_K, swdown_Wm2, latitude_deg, DOY) {
  
  #Convert Tmin from K to C
  tmin <- tmin_K - 273.15
  
  #Convert swdown from W m-2 to MJ/m2/day
  swdown <- swdown_Wm2 * 86400 / 10**6
  
   
  #Convert latitude from degrees to radians
  latitude <- latitude_deg * pi / 180
  
  
  #Stefan-Boltzmann constant [W m-2 K-4] (section 3.2.2)
  sigma <- 5.67*10^(-8)
  
  
  #Actual vapour pressure (eq. 36), in Pa CHECK IF NEED IN mbar !!!!!!!!!
  vap_pressure <- 610.8 * exp((17.27 * tmin) / (237.3 + tmin))
  
  
  #Day angle [radians] (section 3.2.3 fourth equation, also eq. 88)
  Q0 <- (2 * pi * (DOY - 1)) / 365
  
  
  #Solar inclination [radians]
  #(section 3.2.3 third equation)
  delta <- 0.006918 - 0.39912 * cos(Q0) + 0.070257 * sin(Q0) -
           0.006758 * cos(2*Q0) + 0.000907 * sin(2*Q0) -
           0.002697 * cos(3*Q0) + 0.00148 * sin(3*Q0)
  
  
  #Omega: Sunset hour angle [radians] (eq. 87)
  sunset_angle <- acos(-1 * tan(latitude) * tan(delta))
    
  
  #Kd0: expected downwelling shortwave radiation on a cloudless day [MJ m-2 d-1]
  #(section 3.2.3 second equation)
  swdown_clear <- (94.5 * (1 + 0.033 * cos((2*pi*DOY / 365))) / pi) *
                  (sunset_angle * sin(delta) * sin(latitude) + cos(delta) * cos(latitude) * 
                  sin(sunset_angle))
  
  
  #Incoming longwave radiation [W m-2] (section 3.2.3 first equation)
  lwdown <- sigma * tmean_K^4 * (1 - (1- 0.65*(vap_pressure / tmean_K)^0.14) *
                                   (1.35 * (swdown / swdown_clear) - 0.35))
  
  
  return(lwdown)

}



  
  
  
  
  
  
  
  



