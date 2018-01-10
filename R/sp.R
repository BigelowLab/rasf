#' Get one of the common projection strings
#'
#' @export
#' @param name the name of the projection string (just one)
#' @return a projection string
get_proj_string = function(name = c('longlat', "lcc", "all")[1]){
    
    PROJ = c(
        longlat =  "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
        lcc = "+proj=lcc +lat_1=25 +lat_0=25 +lon_0=-95 +k_0=1 +x_0=0 +y_0=0 +a=6367470.21484375 +b=6367470.21484375 +units=km +no_defs")
        
    switch(tolower(name[1]),
        'longlat' = PROJ[['longlat']],
        'lonlat' = PROJ[['longlat']],
        'lcc' = PROJ[['lcc']],
        "all" = PROJ,
        paste('projection name not known:', name[1]))
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