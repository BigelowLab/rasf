

#' Make a raster LUT which is essentially a lookup table that for a given
#' location points to the nearest non-mask value in a mask.
#'
#' Each cell of a lut contains a cell address of the closest
#' cell covered by a non-mask value in the input.  Obviously, cells already
#' covered by non-mask will have it's original cell number.  Cells covered by
#' \code{mask_value} in the input will have a cell number of the closest
#' non-mask cell.
#'
#' Note this assigns the closest cell value which is a 1-d index into 2-d space.
#'
#' @export
#' @param x raster mask
#' @param mask_value numeric the value of masked areas, by default NA
#' @param nonreassigned_value either "cellnumber" which indicates cell number should be used
#'        or some numeric value like NA.
#' @return raster of cell addresses.  Where the input, R, had non-mask values the
#'   cell addresses point to the input cell.  Where the input had mask-valued cells,
#'   the output cell addresses point to the nearest non-mask cells in the input.
make_raster_lut <- function(x = make_dummy_mask(),
                            mask_value = NA_real_,
                            nonreassigned_value = "cellnumber"){

  if (!is_raster_type(x)) stop("Input R must be a raster type")

  # create a matrix with cell numbers (ordered by row top to bottom)
  d <- raster_dim(x)
  #allCell <- matrix(seq_len(d[1]*d[2]), d[1], d[2], byrow = TRUE)
  # create raster of cell numbers and reassign missing values
  m <- seq_len(d[['ncell']])
  if (is_raster(x)){
    R <- raster::raster(matrix(m, d[['ncol']], d[['nrow']], byrow = TRUE), template = x)
  } else {
    R <- terra::rast(nrows = d[['nrow']],
                     ncols = d[['ncol']],
                     nlyrs = 1,
                     xmin = terra::xmin(x),
                     xmax = terra::xmax(x),
                     ymin = terra::ymin(x),
                     ymax = terra::ymax(x),
                     crs = terra::crs(x),
                     vals = m)

  }

  allPts <- as_points(R) %>%
    stats::setNames(c("x", "y", "value"))
  if (is.na(mask_value[1])){
    isna <- is.na(allPts$value)
  } else {
    isna <- dplyr::near(allPts$value, mask_value[1])
    # is any values of x are NA then logical comparisons fail
    # raster should either use NA for mask values OR not have any NAs
    if (any(is.na(isna))){
      warning("one or more input raster values missing - check input and mask_value")
      return(R)
    }
  }

  # if none are NA, then we are done
  if (!any(isna)) return(R)
  R[isna] <- 0  # masked
  if (nonreassigned_value[1] != "cellnumber") R[!isna] <- nonreassigned_value[1]


  ixMask <- allPts$value <= 0

  maskedPts <- allPts %>%
    dplyr::filter(ixMask) %>%
    dplyr::select(.data$x, .data$y) %>%
    as.matrix()
  okPts <- allPts %>%
    dplyr::filter(!ixMask) %>%
    dplyr::select(.data$x, .data$y) %>%
    as.matrix()

   if (is_raster(R)){
     # convert to points and cells
     maskedPts <- raster::rasterToPoints(R, function(x) x <= 0)[,c('x','y')]
     maskedCell <- raster::cellFromXY(R, maskedPts)
     #okPts <- raster::rasterToPoints(R, function(x) x > 0)[,c('x','y')]
   } else {
     maskedPts <- terra::as.points(R)
     maskedCell <- terra::cellFromXY(R, maskedPts)
   }

  # magic
  ix <- RANN::nn2(okPts, maskedPts, k = 1)
  ok <- okPts[ix$nn.idx[,1],]
  # compute the new cell values and assign
  reassignedCell <- if (is_raster(R)){
    raster::cellFromXY(R, ok)
  } else {
    terra::cellFromXY(R, ok)
  }
  R[maskedCell] <- reassignedCell
  R
}


#' Convert raster object to point data
#'
#' @export
#' @param x raster mask as BasicRaster or SpatRaster
#' @param fun function or NULL, to filter points - see \code{\link[raster]{rasterToPoints}}
#' @param ... other arguments passed through
#' @return tibble of
#' \itemize{
#' \item{x ala longitude}
#' \item{y ala latitude}
#' \item{value value at cell}
#' }
as_points <- function(x = make_dummy_mask(), fun = NULL, ...){

  if (!is_raster_type(x)) stop("Input R must be a raster type")
  s <- raster_dim(x)

  if (is_raster(x)){
    v <- raster::values(x)
    r <- raster::xyFromCell(x, seq_len(s[['ncell']])) %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(values = v)
  } else {
    rv <- terra::as.points(x, values = TRUE, na.rm = FALSE)
    v <- terra::values(rv)[,1]
    r <- rv %>%
      terra::geom() %>%
      dplyr::as_tibble() %>%
      dplyr::select(dplyr::all_of(c("x", "y"))) %>%
      dplyr::mutate(value = v)
  }

  if (!is.null(fun)) {
    r <- r %>%
      dplyr::filter(fun(v))
  }

  r
}

#' Determine the closest non-NA pixel given [lon,lat] and a lut raster
#'
#' @export
#' @param x the data frame with lon and lat coordinates
#' @param lut the raster look-up with precomputed closest non-NA cells
#' @return tibble with lon and lat columns
closest_available_cell <- function(
  x = dplyr::tibble(lon = c(0.1, 0.6, 0.9, 0.1),
                    lat = c(0.3, 0.1, 0.8, 0.9)),
  lut = make_raster_lut()){

  if (is_raster(lut)){
    index <- raster::raster(matrix(seq_len(raster::ncell(lut)),
                                   ncol = ncol(lut),
                                   nrow = nrow(lut),
                                   byrow = TRUE),
                            template = lut)
    xy <- x %>% dplyr::select(.data$lon, .data$lat)
    lutCell <- raster::extract(lut, xy)
    indexCell <- raster::extract(index,  xy)
    d <- lutCell != indexCell
    if (any(d)){
      dxy <- raster::xyFromCell(lut, lutCell[d])
      xy$lon[d] <- dxy[,1]
      xy$lat[d] <- dxy[,2]
    }
  } else {
    index <- terra::init(lut, fun = "cell")
    xy <- x %>% dplyr::select(.data$lon, .data$lat)
    lutCell <- terra::extract(lut, xy)[,2]
    indexCell <- terra::extract(index,  xy)[,2]
    d <- lutCell != indexCell
    if (any(d)){
      dxy <- terra::xyFromCell(lut, lutCell[d])
      xy$lon[d] <- dxy[,1]
      xy$lat[d] <- dxy[,2]
    }
  }
  xy
}

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
  if (!is_raster_type(x)) stop("Input x must be a raster type")
  shape <- raster_dim(x)

  if (!is.null(polygon)){
    if (is_raster(x)){
      pcell <- try(raster::extract(x[[1]],
                                   sf::as_Spatial(polygon),
                                   cellnumbers = TRUE)[[1]][,"cell"])
    } else {
      pcell <- try(terra::extract(x[[1]],
                                   terra::vect(sf::as_Spatial(polygon)),
                                   cells = TRUE)[,"cell"])

    }
    if (inherits(pcell, "try-error")){
      stop("unable to extract polygon from raster")
    }
    #nlayers <- seq_len(shape[['nlayer']]-1)
    nlayers <- seq_len(shape[['nlayer']]) - 1
    nindex <- nlayers * shape[['ncell']]
    index <- c(pcell, as.vector(sapply(nindex, function(i) i + pcell)))
    loc <- xyCellLayerFromIndex(index, x)
    values <- extractPts(loc, x)
    loc <- loc %>%
      dplyr::mutate(value = values)
    if (na.rm){
      loc <- loc %>%
        dplyr::filter(!is.na(.data$value))
    }
    if (!is.null(pts)){
      idx <- indexFromPts(pts, x)
      loc <- loc %>%
        dplyr::filter(!(.data$index %in% idx))
    }
    replace <- nrow(loc) < n
    if (replace){
      warning(c(sprintf("number of available cells, %i, less than request, %i ",
                      nrow(loc), n),
                "some points may be duplicated"))

    }
    loc <- loc %>%
      dplyr::sample_n(n, replace = replace)
  } else {
    index <- sample(shape[["nindex"]], n*m, replace = FALSE)
    loc <- xyCellLayerFromIndex(index, x)
    values <- extractPts(loc, x)
    loc <- loc %>%
    	dplyr::mutate(value = values)
    if (na.rm){
      loc <- loc %>%
        dplyr::filter(!is.na(.data$value))
    }
    if (!is.null(pts)){
      if (inherits(pts, 'sf')) pts <- pts %>% sf::st_drop_geometry()
      idx <- indexFromPts(pts, x)
      loc <- loc %>%
        dplyr::filter(!(.data$index %in% idx))
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

  if (inherits(pts, 'sf')) pts <- pts %>% sf::st_drop_geometry()
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
#' @param x raster object
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
  if (!is_raster_type(x)) stop("Input x must be a raster type")
  if (is_raster(x)){
    nrow = as.numeric(raster::nrow(x))
    ncol = as.numeric(raster::ncol(x))
    ncell  = as.numeric(raster::ncell(x))
    nlayer = as.numeric(raster::nlayers(x))
  } else{
    nrow = as.numeric(terra::nrow(x))
    ncol = as.numeric(terra::ncol(x))
    ncell  = as.numeric(terra::ncell(x))
    nlayer = as.numeric(terra::nlyr(x))
  }
  c(ncol = ncol, nrow = nrow, ncell = ncell, nlayer = nlayer, nindex = ncell * nlayer)
}

#' Computer the range of values for raster objects
#'
#' @export
#' @param x a Raster* or SpatRaster class
#' @param na.rm logical, if TRUE remove NAs. For objects inheriting \code{BasicRaster}.
#'   Ignored for class SpatRaster.
#' @param collapse logical, if FALSE then a matrix is returned with two columns
#'    which are vectors of \code{min} and \code{max} values for each layer.
#'    If TRUE then a a two element vector of  \code{min} and \code{max} for all layers
#'    is returned
#' @param ... other arguments for \code{\link[terra]{minmax}}
#' @return either a matrix of min and max values or a two element vector of
#'    min and max
raster_range <- function(x, na.rm = TRUE, collapse = TRUE, ...){
  if (!is_raster_type(x)) stop("Input x must be a raster type")

  if (is_raster(x)){
    mn <- raster::minValue(x)
    mx <- raster::maxValue(x)
    if (collapse[1]){
      r <- c(min(mn, na.rm = na.rm), max(mx, na.rm = na.rm))
    } else {
      r <- cbind(mn, mx)
      colnames(r) <- c("min", "max")
      rownames(r) <- names(x)
    }
  } else {
    r <- terra::minmax(x)
    if (collapse[1]){
      r <- c(min(r[1,], na.rm = na.rm), max(r[2,], na.rm = na.rm))
    } else {
      r <- t(r)
      colnames(r) <- c("min", "max")
      rownames(r) <- names(x)
    }
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
  if (!is_raster_type(x)) stop("Input x must be a raster type")
  shape <- raster_dim(x)
  layer <- ((index-1) %/% shape[['ncell']]) + 1
  cell  <- index - ((layer - 1) * shape[['ncell']])
  col   <- ((cell-1) %% shape[['ncol']])  + 1
  row   <- floor((cell - 1) / shape[['ncol']]) + 1
  if (is_raster(x)){
    xx    <- raster::xFromCol(x, col)
    yy    <- raster::yFromRow(x, row)
  } else {
    xx    <- terra::xFromCol(x, col)
    yy    <- terra::yFromRow(x, row)
  }
    dplyr::tibble(index, cell, col, row, x = xx, y = yy, layer)
}

#' Compute 3d array indices given cell and layer
#'
#' @export
#' @export
#' @param cellLayer 2 element data frame (or tibble) with cell and layer.  If layer is
#'        of type character we covert to integer using \code{names(x)}
#' @param x Raster* of SpatRaster
#' @return cellLayer with index added
indexFromCellLayer <- function(cellLayer, x){
  if (!is_raster_type(x)) stop("Input x must be a raster type")
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
  if (inherits(pts, 'sf')) pts <- pts %>% sf::st_drop_geometry()
  if (!is_raster_type(x)) stop("Input x must be a raster type")
  shape <- raster_dim(x)
  pts <- vetPts(pts)

  if (is_raster(x)){
    col   <- raster::colFromX(x, pts$x)
    row   <- raster::rowFromY(x, pts$y)
  } else {
    col   <- terra::colFromX(x, pts$x)
    row   <- terra::rowFromY(x, pts$y)
  }
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


#' Extract values from a Raster* of SpatRaster object
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

  if (!is_raster_type(x)) stop("Input x must be a raster type")
  shape <- raster_dim(x)


  if (inherits(pts, 'sf')) pts <- pts %>% sf::st_drop_geometry()
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
    if (is_raster(x)){
      xy <- raster::xyFromCell(x, pts[[ic[[1]]]])
    } else {
      xy <- terra::xyFromCell(x, pts[[ic[[1]]]])
    }
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

  if (shape[["nlayer"]] == 1) {
    # single layers are easy
    if(is_raster(x)){
      r <- raster::extract(x, pts[,c(ix,iy)] %>% as.matrix())
    } else {
      r <- terra::extract(x, pts[,c(ix,iy)] %>% as.matrix())
    }
    return(r)
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
                        if (is_raster(x)){
                          r <- raster::extract(x[[j]], pp[[i]][,c(ix, iy), drop = FALSE])
                        } else {
                          r <- terra::extract(x[[j]], pp[[i]][,c(ix, iy), drop = FALSE])[,2]
                        }
                      } else {
                        r <- rep(NA_real_, nrow(pp[[i]]))
                      }
                      r
                      }, simplify = FALSE)
  v       <- unsplit(vv, ff, drop = FALSE)

  unname(v)
}


#' Test if an object is RasterLayer, RasterStack or RasterBrick
#'
#' @export
#' @param x object to test
#' @return logical, TRUE if the object inherits from the specified class
is_raster <- function(x){

  inherits(x, "BasicRaster")
}

#' Test if an object is SpatRaster
#'
#' @export
#' @param x object to test
#' @return logical, TRUE if the object inherits from the specified class
is_terra <- function(x){
  inherits(x, "SpatRaster")
}


#' Test if an object is from \code{raster} or \code{terra} package
#'
#' @export
#' @param x object to test
#' @param klass character a vector of allowed class types.  Use this to narrow
#'   the test, say for a brick by setting klasses to 'RasterBrick'
#' @return logical, TRUE if the object inherits from the specified class
is_raster_type <- function(x,
  klass = c("BasicRaster", "SpatRaster")){
  any(sapply(klass, function(k) {inherits(x, k)}))
}


#' Retrieve metadata about one or more raster files
#'
#' This is a wrapper around \code{\link[rgdal]{GDALinfo}}
#'
#' @export
#' @param x character, one or more filenames (.tiff, .grd, etc.)
#' @return table (tibble) of GDAL metadata
#' \itemize{
#'   \item{rows integer }
#'   \item{columns integer }
#'   \item{bands integer }
#'   \item{ll.x numeric }
#'   \item{ll.y numeric }
#'   \item{res.x numeric }
#'   \item{res.y numeric }
#'   \item{oblique.x numeric }
#'   \item{oblique.y numeric }
#'   \item{crs character }
#' }
raster_fileinfo <- function(x){

  stopifnot(all(file.exists(x)))

  fi <- lapply(x,
    function(f){
      x <- rgdal::GDALinfo(f, returnStats = FALSE, silent = TRUE)
      rbind(c(x[seq_len(length(x))], crs = attr(x, "projection")))
    })
  dplyr::as_tibble(do.call(rbind, fi)) %>%
    dplyr::mutate(rows = as.integer(.data$rows),
                  columns = as.integer(.data$columns),
                  bands = as.integer(.data$bands),
                  ll.x = as.numeric(.data$ll.x),
                  ll.y = as.numeric(.data$ll.y),
                  res.x = as.numeric(.data$res.x),
                  res.y = as.numeric(.data$res.y),
                  oblique.x = as.numeric(.data$oblique.x),
                  oblique.y = as.numeric(.data$oblique.y))
}



#' Crop a raster - a wrapper around \code{raster::crop} or \code{raster::crop}
#'
#' @export
#' @param x RasterLayer or SpatRaster
#' @param ... other arguments for \code{raster::crop} or \code{raster::crop}
#' @return Raster or SpatRaster object
raster_crop <- function(x, ...){

  if (!is_raster_type(x)) stop("Input x must be a raster type")

  r <- NULL

  if (is_terra(x)){
    r <- terra::crop(x, ...)
  } else {
    r <- raster::crop(x, ...)
  }

  r
}


#' Shift a raster - a wrapper around \code{raster::shift} or \code{raster::shift}
#'
#' @export
#' @param x RasterLayer or SpatRaster
#' @param ... other arguments for \code{raster::shift} or \code{raster::shift}
#' @return Raster or SpatRaster object
raster_shift <- function(x, ...){

  if (!is_raster_type(x)) stop("Input x must be a raster type")

  r <- NULL

  if (is_terra(x)){
    r <- terra::shift(x, ...)
  } else {
    r <- raster::shift(x, ...)
  }

  r
}


#' A wrapper around \code{\link[raster]{rotate}} that provides rotation from
#' [-180,180] to [0,360].
#'
#' @export
#' @param x raster, a raster object with global coverage
#' @param inverse logical if FALSE then pass all of the arguments to
#'        \code{\link[raster]{rotate}}.  If FALSE then transform from [-180,180]
#'        to [0,360]
#' @param adjust_origin logical, if TRUE then make minor adjustment to origin.
#'        Ignored if \code{inverse} is FALSE.
#' @param filename character see \code{\link[raster]{merge}}
#' @param ... optional arguments for \code{\link[raster]{merge}}
#' @return raster object
raster_rotate <- function(x,
                          inverse = FALSE,
                          adjust_origin = TRUE,
                          filename = "", ...){
  if (is_raster(x)){
    r <- raster_rotate_raster(x,
                              inverse = inverse,
                              adjust_origin = adjust_origin,
                              filename = filename,
                              ...)
  } else {
    r <- raster_rotate_terra(x,
                              inverse = inverse,
                              adjust_origin = adjust_origin,
                              filename = filename,
                              ...)
  }
  return(r)
}


raster_rotate_terra <- function(x,
                                inverse = inverse,
                                adjust_origin = adjust_origin,
                                filename = filename,
                                ...){
  stop("not implemented for terra (yet)")
}



raster_rotate_raster <- function(x,
                                 inverse = inverse,
                                 adjust_origin = adjust_origin,
                                 filename = filename,
                                 ...) {
  if (inverse){
    # [-180, 180] -> [0,360]
    # modified from original
    # https://github.com/rspatial/raster/blob/master/R/rotate.R
    e <- raster::extent(x)
    xrange <- slot(e, "xmax") - slot(e, "xmin")
    if (xrange < 350 | xrange > 370 | slot(e, "xmin") > -170 | slot(e, "xmax") > 190) {
      if (xrange < 350 | xrange > 370 | slot(e, "xmin") < -190 | slot(e, "xmax") > 190) {
        warning('this does not look like an appropriate object for this function')
      }
    }
    bb <- bb_split(as.vector(e), at = 0)
    west <- raster_crop(x, bb[['bb1']])
    east <- raster_crop(x, bb[['bb2']])
    west2 <- raster_shift(west, dx = 360)
    if (adjust_origin) {

      if (inherits(x, "SpatRaster")){
        stop("not implement for SpatRaster for terra package")

      } else {

        raster::origin(west2) <- raster::origin(west)
      }
    }

    x <- raster::merge(east, west2,
                       filename = filename, ...)

  } else {
    x <- raster::rotate(x, filename = filename, ...)
  }
  x

}
