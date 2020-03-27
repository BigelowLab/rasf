#' Test if WKT is required for specifying coordinate reference
#' 
#' @export
#' @param package character, the name of the package to test
#' @param ... arguments for \code{\link{sf_use_wkt}} or \code{\link{raster_use_wkt}}
#' @return logical TRUE if WKT is required
use_wkt <- function(package = c("sf", "raster")[1], ...){
  switch(tolower(package[1]),
         'raster' = raster_use_wkt(...),
         sf_use_wkt(...))
}
#' Test if we need WKT or proj4string for sf package
#' 
#' @export
#' @param name character, the name of the sf dependency to test
#' @param threshold character, the dependency version that determines if wkt is required
#' @return logical TRUE if WKT is required
sf_use_wkt <- function(name = 'proj.4', threshold = '6'){
  sf::sf_extSoftVersion()[[name]] >= '6'
}


#' Test if we need WKT or proj4string for raster package
#' 
#' @export
#' @param threshold character, the version that determines if wkt is required
#' @return logical TRUE if WKT is required
raster_use_wkt <- function(threshold = '4'){
  packageVersion('raster') >= threshold
}