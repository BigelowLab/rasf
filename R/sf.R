#' Test the geometry for inheritance
#'
#' @export
#' @param x sf object to test
#' @param klass character vector of one or more classes
#' @param geometry_column character the name of the geometry column, by default 'geometry'
#' @return logical
is_sf_geometry <- function(x, klass = "sfc", geometry_column = 'geometry'){
	stopifnot(any(c(inherits(x, "sf"), inherits(x,"sfc"))))
	inherits(x[[geometry_column]], klass)
}

#' Retrieve the class of the geometry
#'
#' @export
#' @param x sf object
#' @param geometry_column character the name of the geometry column, by default 'geometry'
#' @return character
sf_geometry <- function(x, geometry_column = 'geometry'){
	stopifnot(any(c(inherits(x, "sf"), inherits(x,"sfc"))))
	class(x[[geometry_column]])
}


#' Convert sfc_POINTS to a mesh of polygons (delaunay triangulation)
#'
#' @export
#' @param x sf object of type sfc_POINT
#' @param varname character one or more variable names to transfer to the polygon
#' @param fun function to compute value the variables specifed by varname
#' @param ... other arguments for \code{fun}
#' @return sfc_POLYGON table
points_to_mesh <- function(x, varname = NULL, fun = mean, ...){

	stopifnot(is_sf_geometry(x,  "sfc_POINT"))
  xy <- sf::st_coordinates(x)   # n xy locations
  d <- geometry::delaunayn(xy)  # m polygons x 3 nodes
  d <- cbind(d, d[,1])          # wrap the first point
  #iv <- c(1,2,3,1)             # vertex index wraps
  g <- lapply(seq_len(nrow(d)),
  	function(i, xy = NULL, del = NULL){
  			sf::st_polygon(list(xy[del[i,],]))
  	}, xy = xy, del = d)
	mesh <- dplyr::tibble(p1 = d[,1], p2 = d[,2], p3 = d[,3]) %>%
		dplyr::mutate(geometry = g) %>%
    sf::st_sf(sf_column_name = "geometry", crs = sf::st_crs(x)$proj4string)

  if (!is.null(varname)){
  	stopifnot(all(varname %in% colnames(x)))
  	d <- d[,1:3]  # revert to the original m polygons by 3 nodes
  	len <- seq_len(nrow(mesh))  # compute this just once
		for (var in varname){
		  vals <- x[[var]]
			v <- sapply(len,
                  function(i){
                    fun(vals[d[i,]], ...)
                  })
			mesh <- mesh %>%
      	tibble::add_column(!!var := v, .before = "p1")
		}
  }
  mesh
}
