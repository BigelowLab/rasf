#' Convert bounding box [0,360] longitudes to [-180, 180]
#'
#' Bounding boxes are 4 element vectors of [left, right, bottom, top]
#'
#' @export
#' @param x numeric bounding box vector, no check is done for being withing 0,360 range
#' @return numeric bounding box vector
to180BB <- function(x) {x[1:2] <- to180(x[1:2]) ; x}

#' Convert [-180,180] bounding box longitudes to [0,360]
#'
#' Bounding boxes are 4 element vectors of [left, right, bottom, top]
#'
#' @export
#' @param x numeric bounding box vector, no check is done for being withing 0,360 range
#' @return numeric bounding box vector
to360BB <- function(x) {x[1:2] <- to360(x[1:2]) ; x}

#' Convert [0,360] longitudes to [-180, 180]
#'
#' @export
#' @param x numeric vector, no check is done for being withing 0,360 range
#' @return numeric vector
to180 <- function(x) {ix <- x > 180 ; x[ix] <- x[ix]-360; x}

#' Convert [-180,180] longitudes to [0,360]
#'
#' @export
#' @param x numeric vector, no check is done for being withing 0,360 range
#' @return numeric vector
to360 <- function(x) {ix <- x < 0 ; x[ix] <- x[ix]+ 360; x}

#' Get one of the common projection strings
#'
#' @export
#' @param name the name of the projection string (just one)
#' @return a projection string
get_proj_string = function(name = c('longlat', "lcc", "utm-19", "EPSG:3857", "tmerc", "all")[1]){

    PROJ = c(
        "EPSG:3857" = '+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs',
        "longlat" =  "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
        "utm-19" = "+proj=utm +zone=19 +datum=NAD83 +units=m +no_defs",
        "lcc" = "+proj=lcc +lat_1=25 +lat_0=25 +lon_0=-95 +k_0=1 +x_0=0 +y_0=0 +a=6367470.21484375 +b=6367470.21484375 +units=km +no_defs",
        "daymet-lcc" = "+proj=lcc +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +lat_1=25 +ellps=WGS84 +lat_2=45",
        "tmerc" = "+proj=tmerc +lat_0=.... +lon_0=... +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        "merc" = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs")

    switch(tolower(name[1]),
        'longlat' = PROJ[['longlat']],
        'lonlat' = PROJ[['longlat']],
        "utm-19" = PROJ[['utm-19']],
        'lcc' = PROJ[['lcc']],
        "epsg:3857" = PROJ[["EPSG:3857"]],
        "daymet-lcc" = PROJ[["daymet-lcc"]],
        "tmerc" = PROJ[["tmerc"]],
        "merc" = PROJ[['merc']],
        "all" = PROJ,
        paste('projection name not known:', name[1]))
}

#' Get a projection string
#'
#' @seealso \code{link{get_proj_string}}
#' @export
#' @param ... arguments for \code{link{get_proj_string}}
#' @return a projection string
get_crs <- function(...){
    get_proj_string(...)
}



#' Transform locations in a tibble or data frame from one projection to another
#'
#' Elements of the input that are NA are returned as NA
#'
#' @export
#' @param x tibble with columns lon and lat
#' @param from_proj projection of the source coordinates
#' @param to_proj projection of the source coordinates
#' @param from_names string, the names of the x and y coords to project from c('lon', 'lat') default
#' @param to_names string, the names of the x and y coords to project into c('x', 'y') default,
#'        the if from_names is the same as to_names then replacement occurs
#' @return updated tibble
project_tibble <- function(x,
    from_proj = get_proj_string('longlat'),
    to_proj = get_proj_string('lcc'),
    from_names = c('lon', 'lat'),
    to_names = c("x", "y")){

    ll <- x %>% dplyr::select(from_names)
    ix <- apply(ll, 1, function(x) any(is.na(x)) )
    input <- ll %>% filter(!ix)
    sp::coordinates(input) <- from_names
    sp::proj4string(input) <- from_proj
    output <- sp::coordinates(sp::spTransform(input, to_proj))
    ll[!ix, ] <- tibble::as_tibble(output)
    colnames(ll) <- to_names
    if (tibble::has_name(x, to_names[1])) {
       x[[to_names[[1]]]] <- ll[[to_names[1]]]
    } else {
        x <- x %>% tibble::add_column(!!to_names[1] := ll[[to_names[1]]])
    }
    if (tibble::has_name(x, to_names[2])) {
        x[[to_names[[2]]]] <- ll[[to_names[2]]]
    } else {
        x <- x %>% tibble::add_column(!!to_names[2] := ll[[to_names[2]]])
    }
    x
}



#' Convert a 4-element bbox vector to a matrix of two columns (x and y)
#'
#' @export
#' @param x a 4-element numeric vector of left, bottom, right, top coordinates
#' @param close logical, if TRUE then close the polygon such that the first
#'    and last verices are the same
#' @return a matrix of 2 columns and either 5 rows (closed) or 4 rows (open)
bbox_to_matrix <- function(x = c(-72,-63,39,46), close = TRUE){
   if (close) {
     x <-  matrix(c(
            x[1],x[2],x[2],x[1], x[1],
            x[3],x[3],x[4],x[4], x[3]),
           ncol = 2)
   } else {
      x <- matrix(c(
            x[1],x[2],x[2],x[1],
            x[3],x[3],x[4],x[4]),
           ncol = 2)
   }
   x
}

#' Convert a 4-element bbox vector to SpatialPoints
#'
#'
#' @export
#' @param bb a 4-element numeric vector of left, bottom, right, top coordinates
#' @param proj_string a proj4string suitable to pass to \code{sp::CRS()}
#' @return a SpatialPoints object
bbox_to_SpatialPoints<- function(bb = get_bb("gom"),
   proj_string = get_proj_string("longlat")){
   sp::SpatialPoints(bbox_to_matrix(bb), proj4string = sp::CRS(proj_string))
}



#' Convert a 4-element bbox vector to a SpatialPolygons object
#'
#' @export
#' @param bb a 4-element numeric vector of left, bottom, right, top coordinates
#' @param proj_string a proj4string suitable to pass to \code{sp::CRS()}
#' @return a SpatialPolygons object
bbox_to_SpatialPolygons <- function(bb = get_bb("gom"),
   proj_string = get_proj_string("longlat")){
   bb_p <- sp::Polygon(bbox_to_matrix(bb))
   bb_ps <- sp::Polygons(list(bb_p), "bb")
   sp::SpatialPolygons(list(bb_ps), proj4string = sp::CRS(proj_string))
}



#' Convert a 4-element bbox vector to a SpatialPolygonsDataFrame object
#'
#' @export
#' @param bb a 4-element numeric vector of left, bottom, right, top coordinates
#' @param ... further arguments for \code{bbox_to_SpatialPolygons}
#' @return a SpatialPolygons object
bbox_to_SpatialPolygonsDataFrame <- function(bb = get_bb("gom"),...){
   spolys <- bbox_to_SpatialPolygons(bb, ...)
   sp::SpatialPolygonsDataFrame(spolys,
      data = data.frame(ID = names(spolys), row.names = names(spolys)))
}

#' Convert a 4-element bbox vector to a sf::st_bbox object
#'
#' @export
#' @param bb a 4-element numeric vector of left, bottom, right, top coordinates
#' @param proj_string a proj4string suitable to pass to \code{sf::st_crs()}
#' @return a sf:bbox object
bbox_to_sf <- function(bb = get_bb("gom"),
    proj_string = get_proj_string()){

    sf::st_bbox(c(xmin = bb[1], xmax = bb[2], ymin = bb[3], ymax = bb[4]),
                crs = sf::st_crs(proj_string))
}

#' Convert a 4-element bbox vector to a sf::POLYGON simple feature object
#'
#' @export
#' @param bb a 4-element numeric vector of left, bottom, right, top coordinates
#' @param proj_string a proj4string suitable to pass to \code{sf::st_crs()}
#' @return a sf:sf with a single POLYGON geometry object
bbox_to_sfPOLYGON <- function(bb = get_bb("gom"),
    proj_string = get_proj_string()){


    p <- sf::st_polygon(x = list(cbind(bb[c(1,2,2,1,1)], bb[c(3,3,4,4,3)]))) %>%
        sf::st_sfc(crs = proj_string) %>%
        sf::st_sf()

}
