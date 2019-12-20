#' Convert bounding box [0,360] longitudes to [-180, 180]
#'
#' Bounding boxes are 4 element vectors of [left, right, bottom, top]
#'
#' @export
#' @param x numeric bounding box vector, no check is done for being withing 0,360 range
#' @return numeric bounding box vector
to180BB <- function(x) {
	x[1:2] <- to180(x[1:2]) 
	if (identical(x[1], 180)) x[1] <- -180
	x}

#' Convert [-180,180] bounding box longitudes to [0,360]
#'
#' Bounding boxes are 4 element vectors of [left, right, bottom, top]
#'
#' @export
#' @param x numeric bounding box vector, no check is done for being withing 0,360 range
#' @return numeric bounding box vector
to360BB <- function(x) {x[1:2] <- to360(x[1:2]);x}

#' Convert [0,360] longitudes to [-180, 180]
#'
#' @export
#' @param x numeric vector, no check is done for being withing [0, 360] range
#' @return numeric vector
to180 <- function(x) {((x + 180) %% 360) - 180} 

#' Convert [-180,180] longitudes to [0, 360]
#'
#' @export
#' @param x numeric vector, no check is done for being within [0,3 60] range
#' @return numeric vector
to360 <- function(x) {ix <- x < 0 ; x[ix] <- x[ix]+ 360; x}


#' Convert a 4-element bbox vector to a sf bbox object
#'
#' @export
#' @param bb a 4-element numeric vector of left, right, bottom and top coordinates
#' @param crs character, the coordinate reference system
#' @return sf bbox object
bb_to_bbox <- function(bb = c(-72, -63, 39, 46),
                       crs = "+init=epsg:4326"){

  sf::st_bbox(c(xmin = bb[1], xmax = bb[2], ymin = bb[3], ymax = bb[4]),
              crs = crs)
}

#' Convert a 4-element bounding box vector to a polygon simple feature object
#'
#' @export
#' @param bb a 4-element numeric vector of left, right, bottom and top coordinates
#' @param crs character, the coordinate reference system
#' @return a simple feature with a POLYGON
bb_to_polygon <- function(bb = c(-72, -63, 39, 46),
                          crs = "+init=epsg:4326"){


  p <- sf::st_polygon(x = list(cbind(bb[c(1,2,2,1,1)], bb[c(3,3,4,4,3)]))) %>%
    sf::st_sfc(crs= crs) %>%
    sf::st_sf()

}
