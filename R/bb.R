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
	if (identical(x[2], -180)) x[2] <- 180
	x
}

#' Convert [-180,180] bounding box longitudes to [0,360]
#'
#' Bounding boxes are 4 element vectors of [left, right, bottom, top]
#'
#' @export
#' @param x numeric bounding box vector, no check is done for being withing 0,360 range
#' @return numeric bounding box vector
to360BB <- function(x) {
	x[1:2] <- to360(x[1:2])
	if (identical(x[1], 360)) x[1] <- 0   # western edge
	if (identical(x[2], 0)) x[2] <- 360   # eastern edge
	x
}

#' Convert [0,360] longitudes to [-180, 180]
#'
#' @seealso \url{https://gis.stackexchange.com/questions/201789/verifying-formula-that-will-convert-longitude-0-360-to-180-to-180/201793}
#' @export
#' @param x numeric vector, no check is done for being withing [0, 360] range
#' @return numeric vector
to180 <- function(x) {
  x <- ((x + 180) %% 360) - 180
  x
}

#' Convert [-180,180] longitudes to [0, 360]
#'
#' @seealso \url{https://gis.stackexchange.com/questions/201789/verifying-formula-that-will-convert-longitude-0-360-to-180-to-180/201793}
#' @export
#' @param x numeric vector, no check is done for being within [0,3 60] range
#' @return numeric vector
to360 <- function(x) {x %% 360}


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

#' Convert a 4-element bounding box vector to a polygon simple feature object.
#'
#' Alternatively, input a list of bounding boxes, each is transformed into one
#' feature in a multi-feature object.
#'
#' @export
#' @param bb a 4-element numeric vector of left, right, bottom and top
#'        coordinates or a list of the same.
#' @param crs character, the coordinate reference system
#' @return a simple feature with a POLYGON or multiple feature object
bb_to_polygon <- function(bb = c(-72, -63, 39, 46),
                          crs = "+init=epsg:4326"){

  if (is.list(bb)){
    pp <- lapply(seq_len(length(bb)),
                 function(i){
                   x <- bb[[i]]
                   sf::st_polygon(x = list(cbind(x[c(1,2,2,1,1)],
                                                 x[c(3,3,4,4,3)]))) %>%
                     sf::st_sfc(crs = crs) %>%
                     sf::st_sf() %>%
                     dplyr::mutate(ID = i)
                 })
    p <- do.call(rbind, pp)
  } else {
    p <- sf::st_polygon(x = list(cbind(bb[c(1,2,2,1,1)], bb[c(3,3,4,4,3)]))) %>%
      sf::st_sfc(crs= crs) %>%
      sf::st_sf()
  }
  p
}

#' Split a bounding box into two at \code{at}
#'
#' @export
#' @param bb numeric, 4 element bouding box of left, right, bottom and top coordinates
#' @param at numeric, longitude to split around
#' @return list of one or two bounding box vectors
bb_split <- function(bb = c(-170,50,-60,60),
                     at = 0){
  if (bb_straddles(bb, at = at)){
    x <- list(
      bb1 = c(bb[1], at, bb[3:4]),
      bb2 = c(at, bb[2:4])
    )
  } else {
    x <- list( bb1 = bb  )
  }
  x
}

#' Test if a blunding box straddles a longitude
#'
#' @export
#' @param bb numeric, 4 element bouding box
#' @param at numeric, longitude to split around
#' @return logical
bb_straddles <- function(bb = c(-170,50,-60,60),
                         at = 0){
  bb[1] < at && bb[2] > at
}
