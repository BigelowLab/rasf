#' Retrieve a raster of Auckland’s Maunga Whau volcano topography.  Adapted from
#' USGS-R inlmisc package
#'
#' @seealso \url{https://waterdata.usgs.gov/blog/inlmiscmaps/}
#' @seealso \url{https://CRAN.R-project.org/package=inlmisc}
#' @export
#' @param indexed logical, if TRUE then assign 1,2,3,... as cell values
#' @return RasterLayer
volcano_raster <- function(indexed = FALSE){
  m <- t(datasets::volcano)[61:1, ]
  s <- 10
  x <- seq(from = 6478705, length.out = 87, by = s)
  y <- seq(from = 2667405, length.out = 61, by = s)
  r <- raster::raster(m,
                      xmn = min(x) - s/2,
                      xmx = max(x) + s/2,
                      ymn = min(y) - s/2,
                      ymx = max(y) + s/2,
                      crs = "+init=epsg:27200")
  if (indexed){
    s <- raster_dim(r)
    r <- raster::setValues(r, seq_len(s[['ncell']]))
  }
  r
}

#' Retrieve a raster stack of Auckland’s Maunga Whau volcano topography.
#' Layers 2 to nlayers slightly altered from the original as per
#' \code{\link{volcano_raster}}
#'
#' @seealso \code{\link{volcano_raster}}
#' @export
#' @param nlayers numeric, the number of layers for the stack
#' @param indexed logical, if TRUE then assign 1,2,3,... as cell values
#' @return RasterStack
volcano_stack <- function(nlayers = 3,
                          indexed = FALSE){
  r <- volcano_raster(indexed = indexed)
  if (indexed){
      s <- raster_dim(r)
      m <- (seq_len(nlayers) - 1) * s[['ncell']]
      rr <- lapply(seq_len(nlayers), function(i)  r + m[i])
  } else {
    s <- runif(nlayers, min = 0.8, max = 1.2)
    rr <- lapply(seq_len(nlayers), function(i)  r * s[i])
  }
  raster::stack(rr)
}

#' Generate a table of random points in a raster stack (or layer)
#'
#' @export
#' @param x RasterLayer or RasterStack
#' @param n numeric, the numer of points to select
#' @param ... further arguments for \link{randomPts}
#' @return tibble of locations with values
volcano_points <- function(x = volcano_stack(), n = 100, ...){
  randomPts(x, n = 100, ...)
}


#' Generate a polygon for the volcano raster.
#'
#' @export
#' @return simple feature geometry for a POLYGON
#' @examples
#' \dontrun{
#' library(dplyr)
#' x <- volcano_stack(indexed = TRUE)
#' p <- volcano_polygon()
#' pts <- randomPts(x, polygon = p)
#' par(mfrow = c(1,3))
#' for (i in seq_len(3)){
#'   plot(x[[i]], main = paste("Layer", i))
#'   plot(p, add = TRUE)
#'   with(pts %>% dplyr::filter(layer == i), points(x, y))
#' }
#' }
volcano_polygon <- function(){
  x <- 6478700 + c(301, 622, 622, 500, 500, 301, 301)
  y <- 2667400 + c(100, 100, 450, 450, 200, 200, 100)
  sf::st_sfc(sf::st_polygon(list(cbind(x,y))),
               crs = "+init=epsg:27200")
}


