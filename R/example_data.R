#' Create a dummy raster mask for testing \code{make_ratser_lut}
#'
#' @export
#' @param nc integer, number of columns
#' @param nr integer, number of rows
#' @param klass character, one of \code{RasterLayer} of \code{SpatRaster}
#' @return raster mask with NA values assigned to area to be masked
make_dummy_mask <- function(nc = 10, nr = 10,
                            klass = c("RasterLayer", "SpatRaster")[2]){
  m <- matrix(seq_len(nc*nr), ncol = nc, nrow = nr, byrow = TRUE)
  m[lower.tri(m)] <- NA
  if (tolower(klass[1]) == "rasterlayer"){
    r <- raster::raster(m)
  } else {
    r <- terra::rast(t(m))
  }
  r
}

#' Create a dummy raster stack for testing.
#'
#' @export
#' @param nc integer, number of columns
#' @param nr integer, number of rows
#' @param nl integer, number of layers
#' @return RasterStack or multi-layer SpatRaster with some NA cells
make_dummy_stack <- function(nc = 10, nr = 10, nl = 10){
  m <- make_dummy_mask(nc = nc, nr = nr)
  n <- nc * nr
  mm <- lapply(seq_len(nl),
               function(i) {
                 m + ((i-1) * n)
               })
  if (is_raster(m)){
    s <- raster::stack(mm)
    names(s) <- paste("lyr", seq_len(nl), sep =  ".")
  } else {
    s <- Reduce(c, mm)
    names(s) <- paste("lyr", seq_len(nl), sep = ".")
  }
  s
}

#' Make dummy points set for a given raster
#'
#' @export
#' @param R RasterLayer, Stack or Brick or SpatRaster
#' @param N integer, the number of points to generate
#' @param M numeric, a multiplier to start with a selection of random points
#'          from which N are selected.  Ignored if na.rm is \code{FALSE}.
#' @param na.rm logical, if TRUE then avoid NA cells
#' @return tibble of point location info ala \code{\link{extractPts}}
make_dummy_points <- function(R = make_dummy_stack(), N = 10, M = 2, na.rm = TRUE){
  if (!is_raster_type(R)) stop("Input R must be a raster type")
  s <- raster_dim(R)
  if (na.rm){
    index = sample(s["nindex"], N * M, replace = FALSE)
    #pts <- cellLayerFromIndex(R, index)
    pts <- xyCellLayerFromIndex(index, R)
    #v <- layers_extractPoints(R, pts)
    v <- extractPts(pts, R)
    pts <- pts %>%
      dplyr::filter(!is.na(v)) %>%
      dplyr::sample_n(N)
  } else {
    index = sample(s["nindex"], N, replace = FALSE)
    pts <- xyCellLayerFromIndex(index, R)
  }
  pts
}



#' Retrieve a raster of Auckland’s Maunga Whau volcano topography.  Adapted from
#' USGS-R inlmisc package.
#'
#' @seealso \url{https://waterdata.usgs.gov/blog/inlmiscmaps/}
#' @seealso \url{https://CRAN.R-project.org/package=inlmisc}
#' @export
#' @param indexed logical, if TRUE then assign 1,2,3,... as cell values
#' @param klass character, one of \code{RasterLayer} of \code{SpatRaster}
#' @return RasterLayer or SpatRaster
volcano_raster <- function(indexed = FALSE,
                           klass = c("RasterLayer", "SpatRaster")[2]){

  m <- t(datasets::volcano)[61:1, ]
  s <- 10
  x <- seq(from = 6478705, length.out = 87, by = s)
  y <- seq(from = 2667405, length.out = 61, by = s)
  
  if (tolower(klass[1]) == 'rasterlayer'){
    crs <- ifelse(use_wkt("raster"), "epsg:27200", "+init=epsg:27200")
    r <- raster::raster(m,
                        xmn = min(x) - s/2,
                        xmx = max(x) + s/2,
                        ymn = min(y) - s/2,
                        ymx = max(y) + s/2,
                        crs = crs)
    if (indexed){
      s <- raster_dim(r)
      r <- raster::setValues(r, seq_len(s[['ncell']]))
    }
  } else {
    crs <- ifelse(use_wkt("terra"), "epsg:27200", "+init=epsg:27200")
    r <- terra::rast(nrows = nrow(m),
                     ncols = ncol(m),
                     xmin = min(x) - s/2,
                     xmax = max(x) + s/2,
                     ymin = min(y) - s/2,
                     ymax = max(y) + s/2,
                     crs = crs,
                     vals = as.vector(t(m)))
    if (indexed){
      s <- raster_dim(r)
      r <- terra::setValues(r, seq_len(s[['ncell']]))
    }                 
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
#' @param ... further arguments for \code{\link{volcano_raster}}
#' @return RasterStack or SpatRaster
volcano_stack <- function(nlayers = 3,
                          indexed = FALSE, 
                          ...){
  r <- volcano_raster(indexed = indexed, ...)
  if (inherits(r, "SpatRaster")){
    
    if (indexed){
      s <- raster_dim(r)
      m <- (seq_len(nlayers) - 1) * s[['ncell']]
      rr <- lapply(seq_len(nlayers), function(i)  r + m[i])
    } else {
      s <- runif(nlayers, min = 0.8, max = 1.2)
      rr <- lapply(seq_len(nlayers), function(i)  r * s[i])
    }
    r <- Reduce(c, rr)
    
  } else {
    if (indexed){
        s <- raster_dim(r)
        m <- (seq_len(nlayers) - 1) * s[['ncell']]
        rr <- lapply(seq_len(nlayers), function(i)  r + m[i])
    } else {
      s <- runif(nlayers, min = 0.8, max = 1.2)
      rr <- lapply(seq_len(nlayers), function(i)  r * s[i])
    }
    r <- raster::stack(rr)
  }
  r
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
  crs <- ifelse(use_wkt("sf"), "epsg:27200", "+init=epsg:27200")
  x <- 6478700 + c(301, 622, 622, 500, 500, 301, 301)
  y <- 2667400 + c(100, 100, 450, 450, 200, 200, 100)
  sf::st_sfc(sf::st_polygon(list(cbind(x,y))),
               crs = crs)
}

#' Generate a dataset of points for use with \code{\link{st_hexbin}}
#'
#' @export
#' @param n the number of points to generate
#' @param bb 4 element numeric of bounding box [left, right, bottom, top]
hexbin_points <- function( n = 1000, bb = c(-72, -63, 39, 46)){
  bb = c(-72, -63, 39, 46)
  xr <- range(bb[1:2])
  yr <- range(bb[3:4])
   dplyr::tibble(
     x = runif(n, min = bb[1], max = bb[2]),
     y = runif(n, min = bb[3], max = bb[4]),
     z = runif(n, min = 0, max = 1))
}

