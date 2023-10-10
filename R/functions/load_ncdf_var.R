#Opens netcdf, retrieves variable and closes file connection
load_ncdf_var   <- function(file, var)
{
  out_ncdf <- nc_open(file)
  out_data <- ncvar_get(out_ncdf, varid=var)
  nc_close(out_ncdf)
  
  return(out_data)
}
