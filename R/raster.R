#
# if polygon
#    pool =  all points for all layers in polygon
#    if na.rm
#       pool = remove missing
#    if pts
#       pool = remove coincident
#    pool = switch
#       sample if enough points left
#       all remaing with warning
#
# else
#    pool = random sample of n*m
#    if na.rm
#       pool = remove missing
#    if pts
#       pool = remove coincident

#' Generate random points within a raster object.
#'
#' Providing a polygon to limit the search area can speed up the search
#' if the polygon is small relative to the area of a layer.
#'
#' @export
#' @param x Raster* layer, brick or stack
#' @param n numeric, the number of points to return
#' @param m multiplier to use when seeking to avoid missing values.
#'        Ignored if polygon is provided.
#' @param na.rm logical, if TRUE then avoid cells with missing values.
#' @param pts a table of points to avoid. If na.rm is TRUE that is handled first,
#'        then point avoidance is handled.  Ignored if NULL.
#' @param polygon polygon (sfc or SpatialPolygon) that describes the region to select
#'        points from.  Ignored if NULL.
#' @return a table of locations with values.  Note that it is possible
#'         to filter the available pool of candidate cells to something
#'         less than the requested sample, n.  In such cases some
#'         samples may be replicates.
randomPts <- function(x,
                      n = 100,
                      m = 2,
                      na.rm = FALSE,
                      pts = NULL,
                      polygon = NULL){

  shape <- raster_dim(x)

  if (!is.null(polygon)){
    pcell <- try(raster::extract(x[[1]],
                                 sf::as_Spatial(polygon),
                                 cellnumbers = TRUE)[[1]][,"cell"])
    if (inherits(pcell, "try-error")){
      stop("unable to extract polygon from raster")
    }
    nlayers <- seq_len(shape[['nlayer']]-1)
    nindex <- nlayers * shape[['ncell']]
    index <- c(pcell, as.vector(sapply(nindex, function(i) i + pcell)))
    loc <- xyCellLayerFromIndex(index, x)
    loc$value <- extractPts(loc, x)
    if (na.rm){
      loc <- loc %>%
        dplyr::filter(!is.na(.data$value))
    }
    if (!is.null(pts)){
      idx <- indexFromPts(pts, x)
      loc <- loc %>%
        dplyr::filter(!(idx %in% .data$index))
    }
    replace <- nrow(loc) < n
    if (replace){
      warning(c(sprintf("number of available cells, %i, less than request, %i",
                      nrow(loc), n),
                "some points may be duplicated"))

    }
    loc <- loc %>%
      dplyr::sample_n(n, replace = replace)
  } else {
    index <- sample(shape[["nindex"]], n*m, replace = FALSE)
    loc <- xyCellLayerFromIndex(index, x)
    loc$value <- extractPts(loc, x)
    if (na.rm){
      loc <- loc %>%
        dplyr::filter(!is.na(.data$value))
    }
    if (!is.null(pts)){
      idx <- indexFromPts(pts, x)
      loc <- loc %>%
        dplyr::filter(!(idx %in% .data$index))
    }
    replace <- nrow(loc) < n
    if (replace){
      warning(c(sprintf(paste("number of pooled cells, %i, less than request %i",
                            "- try increasing value of m, %i?"),
                      nrow(loc), n, m),
                "some points may be duplicated"))
    }
    loc <- loc %>%
      dplyr::sample_n(n, replace = replace)
  }
  loc
}

#' Vet points to ensure uniform input
#'
#' The following are synonyms and definitions.  Note that layer is not required,
#' if not present then it is assigned the value 1 for all points.
#' \itemize{
#'    \item{x: x, lon}
#'    \item{y: y, lat}
#'    \item{layer: layer, z, name}
#'}
#' @export
#' @param pts a (column-named) matrix, data.frame or tibble
#' @return a tibble of x, y and layer
vetPts <- function(pts){

  if (!(is.data.frame(pts) || is.matrix(pts)))
    stop("pts must be data.frame or matrix")
  if (inherits(pts, "matrix") || inherits(pts, 'tbl_df'))
    pts <- as.data.frame(pts, stringsAsFactors = FALSE)

  nm <- names(pts)
  iz <- which(nm %in% c("layer", "z", "name"))[1]
  # it is not an error to pass a single layer
  #if (length(iz) == 0) stop("pts must have 'layer', 'z' or 'name' column")

  ix <- which(nm %in% c("x", "lon"))[1]
  if (length(ix) == 0) stop("pts must have 'x' or 'lon' column")

  iy <- which(nm %in% c("y", "lat"))[1]
  if (length(iy) == 0) stop("pts must have 'y' or 'lat' column")

  r <- dplyr::tibble(x = pts[,ix[1]], y = pts[,iy[1]])
  if (length(iz > 0)){
    r <- r %>%
      dplyr::mutate(layer = pts[,iz[1]])
  } else {
    r <- r %>%
      dplyr::mutate(layer = 1)
  }
  r
}

#' Assemble a vector of raster dimension
#'
#' @export
#' @param x Raster object
#' @return numeric vector of
#' \itemize{
#'   \item{ncol the number of columns}
#'   \item{nrow the number of rows}
#'   \item{ncell the number of cells (pixels) in each layer (ncol * nrow)}
#'   \item{nlayer the number of layers (bands)}
#'   \item{nindex the number of indices (ncell*nlayer)}
#' }
#' @examples
#' \dontrun{
#' slogo <- stack(system.file("external/rlogo.grd", package="raster"))
#' raster_dim(slog)
#' }
raster_dim <- function(x){
  nrow = as.numeric(raster::nrow(x))
  ncol = as.numeric(raster::ncol(x))
  ncell  = as.numeric(raster::ncell(x))
  nlayer = as.numeric(raster::nlayers(x))
  c(ncol = ncol, nrow = nrow, ncell = ncell, nlayer = nlayer, nindex = ncell * nlayer)
}


#' Computer the range of values for Raster* objects
#'
#' @export
#' @param x a Raster* class
#' @param na.rm logical, if TRUE remove NAs
#' @param collapse logical, if FALSE then a matrix is returned with two columns
#'    which are vectors of \code{min} and \code{max} values for each layer.
#'    If TRUE then a a two element vector of  \code{min} and \code{max} for all layers
#'    is returned
#' @return either a matrix of min and max values or a two element vector of
#'    min and max
raster_range <- function(x, na.rm = TRUE, collapse = TRUE){

  mn <- raster::minValue(x)
  mx <- raster::maxValue(x)

  if (collapse[1]){
    r <- c(min(mn, na.rm = na.rm), max(mx, na.rm = na.rm))
  } else {
    r <- cbind(mn, mx)
    colnames(r) <- c("min", "max")
    rownames(r) <- names(x)
  }
  r
}



#' Convert a index (1, 2, 3,... ncells*nlayers) into a multilayer raster to row, col, cell and
#'    layer coordinates
#'
#' @export
#' @param index vector of multilayer index coordinates
#' @param x Raster* layer, brick or stack
#' @return a tibble of index, cell, col, row, x, y, and layer
xyCellLayerFromIndex <- function(index, x){
  shape <- raster_dim(x)
  layer <- ((index-1) %/% shape[['ncell']]) + 1
  cell  <- index - ((layer - 1) * shape[['ncell']])
  col   <- ((cell-1) %% shape[['ncol']])  + 1
  row   <- floor((cell - 1) / shape[['ncol']]) + 1
  xx    <- raster::xFromCol(x, col)
  yy    <- raster::yFromRow(x, row)
  dplyr::tibble(index, cell, col, row, x = xx, y = yy, layer)
}

#' Compute array indices given cell and layer
#'
#' @export
#' @export
#' @param cellLayer 2 element data frame (or tibble) with cell and layer.  If layer is
#'        of type character we covert to integer using \code{names(x)}
#' @param x Raster* layer, brick or stack
#' @return cellLayer with index added
indexFromCellLayer <- function(cellLayer, x){

  if (!all(c("cell", "layer") %in% names(cellLayer))){
    stop("cell layer must contain both cell and layer")
  }
  shape <- raster_dim(x)
  cell <- cellLayer$cell
  layer <- cellLayer$layer
  layerIsNames <- inherits(layer, "character")
  if (layerIsNames){
    nm <- names(x)
    if (is.null(nm)) stop("if layer is of type character then x must have layer names")
    ilayer <- match(layer, nm)
    if (any(is.na(ilayer))) stop("one or more layers not found:", paste(layer(is.na(ilayer)), collapse = " "))
    layer <- ilayer
  }

  index <- cell + (layer-1) * shape[['ncell']]

  dplyr::mutate(cellLayer, index)
}


#' Compute array indices from points - points are in the raster's CRS system
#'
#' @export
#' @param pts location info for points - see \code{vetPts()}
#' @param x layer, brick or stack
#' @return vector of array indices as if the raster were a 3d array
indexFromPts <- function(pts, x ){
  if (!inherits(x, 'BasicRaster')) stop("Input x must be a Raster* class")
  shape <- raster_dim(x)
  pts <- vetPts(pts)
  col   <- raster::colFromX(x, pts$x)
  row   <- raster::rowFromY(x, pts$y)
  z   <- pts$layer
  if (is.character(z)) {
    nm <- names(x)
    if (is.null(nm)) stop("if layer is character then x must have layer names")
    z <- match(z, names(x))
  }
  #        row + number of rows * rowsize + number of layers * layersize
  #index <- x + (y-1)*nx  + (z-1)*nc
  index <- col + (row-1)*shape[['ncol']] + (z-1)*shape[['ncell']]
  return(index)
}


#' Extract values from a Raster* object
#'
#' @export
#' @param pts location info for points to be extracted. Must be a data frame
#'  or matrix with either set of the following columns. Note that layer (or z)
#'  may be either a layer index number or layer names. If 'cell' is present in the
#'  input then [x or lon] and [y or lat] are ignored.
#' \itemize{
#'      \item{'cell' and 'layer' (or 'z')}
#'      \item{'row', 'col' and 'layer' (or 'z')}
#'      \item{'lon', 'lat' and 'layer' (or 'z')}
#'      \item{'x', 'y' and 'layer' (or 'z')}
#'  }
#' @param x a single or multilayer Raster* object.  If \code{x} has only one layer
#'        then pts is passed to \code{raster::extract()} immediately.
#' @return a vector of values
extractPts <- function(pts, x){

  if (!inherits(x, 'BasicRaster')) stop("Input x must be a Raster* class")
  shape <- raster_dim(x)
  if (shape[["nlayer"]] == 1) return(raster::extract(x, pts))

  if (!inherits(pts, "tbl_df")) {
    pts <- try(dplyr::as_tibble(pts))
    if (inherits(pts, 'try-error')){
      stop("unable to cast pts to tibble")
    }
  }
  pts <- vetPts(pts)
  nm <- names(pts)
  iz <- which(nm %in% c("layer", "z"))[1]
  if (length(iz) == 0) stop("pts must have 'layer' or 'z' column")
  layer <- pts[,iz]

  ic <- which(nm %in% "cell")
  if (length(ic) > 0){
    # if the user provides 'cell' then use that first
    xy <- raster::xyFromCell(x, pts[[ic[[1]]]])
    pts <- dplyr::tibble(
      x = xy[,1],
      y = xy[,2],
      layer = pts[,iz])
    ix <- 1
    iy <- 2
    iz <- 3
  } else {
    # otherwise we must have [x or lon] and [y or lat]
    ix <- which(nm %in% c("x", "lon"))[1]
    if (length(ix) == 0) stop("pts must have 'x' or 'lon' column")

    iy <- which(nm %in% c("y", "lat"))[1]
    if (length(iy) == 0) stop("pts must have 'y' or 'lat' column")
  }
  if (inherits(pts$layer, 'character')){
    nm <- names(x)
    if (is.null(nm)){
      stop("if pts layer is character then raster must have layer names ")
    }
    ilayer <- match(pts$layer, nm)
    if (any(is.na(ilayer))){
      warning("one or more layers in pts not found in raster:",
              paste(pts$layer[is.na(ilayer)], collapse = ", "))
    }
    pts$layer <- ilayer
  }
  ff      <- pts$layer
  pp      <- split(pts, ff)
  layers  <- seq.int(shape[['nlayer']])
  vv      <- sapply(names(pp),
                    function(i){
                      j <- as.integer(i)
                      if (j %in% layers){
                        r <- raster::extract(x[[j]], pp[[j]][,c(ix, iy), drop = FALSE])
                      } else {
                        r <- rep(NA_real_, nrow(pp[[j]]))
                      }
                      r
                      })
  v       <- unsplit(vv, ff, drop = FALSE)

  unname(v)
}


