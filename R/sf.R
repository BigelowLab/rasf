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
  iv <- c(1,2,3,1)  # vertex index wraps
  g <- lapply(seq_len(nrow(d)),
  	function(i, xy = NULL, del = NULL, rowindex = NULL){
  			sf::st_polygon(list(xy[del[i, rowindex],]))
  	}, xy = xy, del = d, rowindex = iv)
	mesh <- dplyr::tibble(p1 = d[,1], p2 = d[,2], p3 = d[,3]) %>%
		dplyr::mutate( geometry = g) %>%
    sf::st_sf(sf_column_name = "geometry", crs = sf::st_crs(x)$proj4string)

  if (!is.null(varname)){
  	stopifnot(all(varname %in% colnames(x)))
		for (var in varname){
			mesh <- mesh %>%
      	tibble::add_column(!!var :=
                    sapply(seq_len(nrow(mesh)),
                           function(i){
                             fun((x %>% dplyr::slice(d[i,]))[[var]], ...)
                           }), .before = "p1")
		}
  }
  mesh
}
