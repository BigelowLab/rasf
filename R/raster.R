#' Determine the closest ocean pixel given [lon,lat] and a lut raster
#'
#' @export
#' @param x the data frame with lon and lat coordinates
#' @param lut the raster look-up with precomputed closest ocean cells
#' @return tibble with lon and lat columns
closest_ocean_cell <- function(x, lut = "/mnt/ecocast/coredata/obpg/world/lut_4km.grd"){
    reassignCell <- raster::extract(lut, x)
    xy <- raster::xyFromCell(lut, reassignCell)
    dplyr::tibble(lon = xy[,1], lat = xy[,2])
}

#' Read the header portion of a raster file
#'
#' @export
#' @param filename on or more filenames
#' @param form return form can be 'list' or 'tibble' (default)
#' @return one or more grd lists or the same as a tibble
#' \itemize{
#'      \item{dim vector of nrows, ncols, ncells, nlayers}
#'      \item{extent vector of xmin, xmax, ymin, ymax}
#'      \item{projection  string}
#'      \item{filename string}
#'  }
read_grd <- function(
    filename = file.path('/mnt/ecocast/coredata/obpg/gom/DAY/chlor_a',
                         'A2014365.L3m_DAY_CHL_chlor_a_4km.grd'),
    form = c("list", "tibble")[2]){

    if (length(filename) > 1) {
        x = lapply(filename, read_grd, form = form)
        if (tolower(form[1]) == 'tibble') x = dplyr::bind_rows(x)
        return(x)
    }

    x <- raster::readIniFile(filename, aslist = TRUE)
    g <- x$georeference
    nc <- as.numeric(g$ncols)
    nr <- as.numeric(g$nrows)
    nl <- as.numeric(x$data$nbands)
    if (tolower(form[1]) == 'list'){
        r = list(
            dim    = c(nr, nr, nc*nr, nl),
            extent = as.numeric(c(g$xmin, g$xmax, g$ymin, g$ymax)),
            projection = g$projection,
            filename = filename)
    } else {
        r = dplyr::tibble(
            nrow        = nr,
            ncol        = nc,
            ncell       = nr*nc,
            nlayer      = nl,
            xmin        = g$xmin,
            xmax        = g$xmax,
            ymin        = g$ymin,
            ymax        = g$ymax,
            projection  = g$projection,
            filename    = filename)
    }
    r
}

#' Extract the dimensions from one or more raster header lists
#'
#' @export
#' @param x one or more grd lists
#' @return a matrix of length(x) rows with columns nrows, ncols, ncells, nlayers
grd_dim <- function(x){
    x <- t(sapply(x, '[[', 'dim'))
    colnames(x) = c("nrows", "ncols", "ncells", "nlayers")
    x
}


#' Vet ponts to ensure uniform input
#'
#' The following are synonyms and definitions
#' \itemize{
#'    \item{x: x, lon}
#'    \item{y: y, lat}
#'    \item{layer: layer, z, name}
#'}
#' @export
#' @param pts a matrix, data.frame or tibble
#' @return a tibble of x, y and layer
vet_pts <- function(pts){

    if (!(is.data.frame(pts) || is.matrix(pts)))
        stop("pts must be data.frame or matrix")
    if (inherits(pts, "matrix") || inherits(pts, 'tbl_df'))
        pts <- as.data.frame(pts, stringsAsFactors = FALSE)

    nm <- names(pts)
    iz <- which(nm %in% c("layer", "z", "name"))[1]
    if (length(iz) == 0) stop("pts must have 'layer', 'z' or 'name' column")

    ix <- which(nm %in% c("x", "lon"))[1]
    if (length(ix) == 0) stop("pts must have 'x' or 'lon' column")

    iy <- which(nm %in% c("y", "lat"))[1]
    if (length(iy) == 0) stop("pts must have 'y' or 'lat' column")

    dplyr::tibble(x = pts[,ix], y = pts[,iy], layer = pts[,iz])
}


#' Assemble a vector of raster dimension
#'
#' @export
#' @param x Raster object
#' @return numeric vector of nc = ncols, nr = nrows, n = ncells (2d),
#'     nl = nlayers, N = total number of cells (3d)
raster_shape <- function(x){
    nr = as.numeric(raster::nrow(x))
    nc = as.numeric(raster::ncol(x))
    n  = as.numeric(raster::ncell(x))
    nl = as.numeric(raster::nlayers(x))
    c(nc = nc, nr = nr, n = n, nl = nl, N = nc*nr*nl)
}
#' Convert a index (1... ncells*nlayers) into a multilayer raster to cell and
#'    layer coordinates
#'
#' @export
#' @param index vector of multilayer index coordinates
#' @param R brick or stack
#' @param shape the shape vector to operate with, by default computed on input R
#' @return a tibble of index, cell, col, row, x, y, and layer
cellLayerFromIndex <- function(R, index, shape = raster_shape(R)){

    # nr = as.numeric(raster::nrow(R))
    # nc = as.numeric(raster::ncol(R))
    # n  = as.numeric(raster::ncell(R))
    # nl = as.numeric(raster::nlayers(R))
    # N  = nc*nr*nl
    shape = raster_shape(R)
    layer  = ((index-1) %/% shape[['n']]) + 1
    cell = index - ((layer - 1) * shape[['n']])
    ix = ((cell-1) %% shape[['nc']])  + 1
    iy = floor((cell - 1) / shape[['nc']]) + 1
    x = raster::xFromCol(R, ix)
    y = raster::yFromRow(R, iy)
    dplyr::tibble(index, cell, col = ix, row = iy, x, y, layer)
}

#' Compute array indices from points
#'
#' @export
#' @param R brick or stack
#' @param pts location info for points - see \code{vet_pts()}
#' @return vector of array indices as if the raster were a 3d array
indexFromPts <- function(R, pts){
    if (!inherits(R, 'BasicRaster')) stop("Input R must be a Raster* class")
    #nc <- as.integer(raster::ncell(R))
    #nl <- as.integer(raster::nlayers(R))
    #ny <- raster::nrow(R)
    #nx <- raster::ncol(R)
    s <- raster_shape(R)

    pts <- vet_pts(pts)
    # rows
    x   <- raster::colFromX(R, pts$x)
    # cols
    y   <- raster::rowFromY(R, pts$y)
    z   <- pts$layer
    if (is.character(z)) z <- match(z, names(R))
    # row + number of rows * rowsize + number of layers * layersize
    #index <- x + (y-1)*nx  + (z-1)*nc
    index <- x + (y-1)*s[['nc']] + (z-1)*s[['n']]
    return(index)
}

#' Converts point location to cell number.
#'
#' This function is essentially a wrapper around the \code{raster::cellFrom*()} functions.
#'
#' @export
#' @param R a multilayer Raster* object
#' @param pts location info for points must have at least one of the following
#'  combinations of columns
#'  \itemize{
#'      \item{'cell' if present this is returned}
#'      \item{'row' and 'col'}
#'      \item{'lon' and 'lat'}
#'      \item{'x' and 'y' assumed to be the same as 'lon' and 'lat'}
#'  }
#' @return vector of cell numbers
cellFromPts <- function(R, pts){

    if (!inherits(R, 'BasicRaster')) stop("Input R must be a Raster* class")
    if (!is.data.frame(pts) && !is.matrix(pts))
        stop("Input pts must be data.frame or matrix")
    if (inherits(pts, 'tbl_df')) pts <- as.data.frame(pts, stringsAsFactors = FALSE)
    # setup
    nm <- colnames(pts)
    cell <- NULL

    # easy
    if ('cell' %in% nm) cell <- pts[['cell']]

    # row/col is the next easiest
    if (!is.null(cell)){
        if (all(c('row', 'col') %in% nm)){
            cell <- raster::cellFromRowCol(R, pts[['row']], pts[['col']])
        }
    }

    # real coordinates provided
    if (is.null(cell)){
        ix <- which(nm %in% c("x", "lon"))[1]
        if (length(ix)> 0){
            iy <- which(nm %in% c("y", "lat"))[1]
            if (length(iy) > 0) cell <- raster::cellFromXY(R, pts[,c(ix,iy)])
        }
    }

    if (is.null(cell)) warning("unable to compute cell - returning NULL")
    return(cell)
}


#' Converts point location to cell number and layer number.
#'
#' This function is essentially a wrapper around the \code{raster::cellFrom*()} functions.
#'
#' @export
#' @param R a multilayer Raster* object
#' @param pts location info for points must have at least one of the following
#'  combinations of columns
#'  \itemize{
#'      \item{'cell' if present this is returned}
#'      \item{'row' and 'col'}
#'      \item{'lon' and 'lat'}
#'      \item{'x' and 'y' assumed to be the same as 'lon' and 'lat'}
#'  }
#' @return vector of cell numbers
cellLayerFromPts <- function(R, pts){

    if (!inherits(R, 'BasicRaster')) stop("Input R must be a Raster* class")
    if (!is.data.frame(pts) && !is.matrix(pts))
        stop("Input pts must be data.frame or matrix")
    if (inherits(pts, 'tbl_df')) pts <- as.data.frame(pts, stringsAsFactors = FALSE)
    # setup
    nm <- colnames(pts)
    cell <- NULL

    # easy
    if ('cell' %in% nm) cell <- pts[['cell']]

    # row/col is the next easiest
    if (!is.null(cell)){
        if (all(c('row', 'col') %in% nm)){
            cell <- raster::cellFromRowCol(R, pts[['row']], pts[['col']])
        }
    }

    # real coordinates provided
    if (is.null(cell)){
        ix <- which(nm %in% c("x", "lon"))[1]
        if (length(ix)> 0){
            iy <- which(nm %in% c("y", "lat"))[1]
            if (length(iy) > 0) cell <- raster::cellFromXY(R, pts[,c(ix,iy)])
        }
    }

    if (is.null(cell)) warning("unable to compute cell - returning NULL")
    return(cell)
}

#' Extract values from a multilayer Raster*
#'
#' @export
#' @param R a multilayer Raster* object
#' @param pts location info for points to be extracted. Must be a data frame
#'  or matrix with either set of the following columns. Note that layer (or z)
#'  may be either a layer index number or layer names.
#' \itemize{
#'      \item{'cell' and 'layer' (or 'z')}
#'      \item{'row', 'col' and 'layer' (or 'z')}
#'      \item{'lon', 'lat' and 'layer' (or 'z')}
#'      \item{'x', 'y' and 'layer' (or 'z')}
#'  }
#' @return a vector of values
layers_extractPoints <- function(R, pts){

    if (!inherits(R, 'BasicRaster')) stop("Input R must be a Raster* class")
    nc <- as.numeric(raster::ncell(R))
    nl <- as.numeric(raster::nlayers(R))
    ny <- as.numeric(raster::nrow(R))
    nx <- as.numeric(raster::ncol(R))

    if (!(is.data.frame(pts) || is.matrix(pts)))
        stop("pts must be data.frame or matrix")
    if (inherits(pts, "matrix") || inherits(pts, 'tbl_df'))
        pts <- as.data.frame(pts, stringsAsFactors = FALSE)

    nm <- names(pts)
    iz <- which(nm %in% c("layer", "z"))[1]
    if (length(iz) == 0) stop("pts must have 'layer' or 'z' column")
    layer <- pts[,iz]
    #if (is.character(layer)) layer <- match(layer, names(R))

    ix <- which(nm %in% c("x", "lon"))[1]
    if (length(ix) == 0) stop("pts must have 'x' or 'lon' column")

    iy <- which(nm %in% c("y", "lat"))[1]
    if (length(iy) == 0) stop("pts must have 'y' or 'lat' column")

    #flayer  <- factor(layer)
    pp      <- split(pts, layer)
    nmR     <- names(R)
    vv      <- sapply(names(pp),
        function(i){
            if (i %in% nmR)
                raster::extract(R[[i]], pp[[i]][,c(ix, iy), drop = FALSE])
            else
                rep(NA_real_, nrow(pp[[i]]))
        })
    v       <- unsplit(vv, layer, drop = FALSE)

    #cell <- cellFromPts(R, pts)

    #index <- cell + (layer-1) * nc

    #raster::getValues(R)[index]
    unname(v)
}



#' Select N random points from a mulitlayer Raster* object
#'
#' This function does not guarantee that returned locations do not coincide with
#' missing points.  The first layer in R is used with \code{dismo::randomPoints()}
#' as if the first layer is a mask so select M*N non-NA points.
#' Those selected points are then randomly distributed across the multiple layers
#' in R - which gives rise to the possibility that NA points may be selected.
#'
#' @export
#' @param R a multilayer Raster* object
#' @param layers a vector of layer IDs.  If not NULL then \code{faux_randomPoints()}
#'  is called instead.
#' @param pts xyz values for presence points [lon, lat, layer] - these
#'  locations are avoided when sampling. Ignored if NULL.   This must have
#'  either set of the following columns.  Note that layer (or z)
#'  may be either a layer index number or layer names.
#' \itemize{
#'      \item{'cell' and 'layer' (or 'z')}
#'      \item{'row', 'col' and 'layer' (or 'z')}
#'      \item{'lon', 'lat' and 'layer' (or 'z')}
#'      \item{'x', 'y' and 'layer' (or 'z')}
#'  }
#' @param N the number of points to select
#' @param M numeric, a multiplier to start with a selection of random points
#'          from which N are selected.
#' @return data.frame of sampled points with the following
#'  \itemize{
#'      \item{lon}
#'      \item{lat}
#'      \item{col}
#'      \item{row}
#'      \item{cell 2D matrix index}
#'      \item{index 3d array index}
#'      \item{layer}
#'      \item{value}
#'  }
layers_randomPoints <- function(R,
    layers = NULL,
    pts = NULL,
    N = 1000,
    M = 5){


    if (!is.null(layers)) {
        x = faux_randomPoints(R, layers = layers, pts = pts,N = N, M = M)
        return(x)
    }

    s   = raster_shape(R)
    # allow dismo::randomPoints to select from the 2d plane of the first layer
    xy  = dismo::randomPoints(R, n = N*M)
    # generate random layers
    z   = sample(seq_len(s[['nl']]), size = N*M, replace = TRUE)
    layers = names(R)[z]
    xyz = dplyr::tibble(x = xy[,'x'], y = xy[,'y'], layer = layers)
    # get the 3d index lof each
    index = indexFromPts(R, xyz)
    if (!is.null(pts)){
        pindex = indexFromPts(R, pts)
        pindex = pindex[!is.na(pindex)]
        index  = index[!(index %in% pindex)]
    }
    index = sample(index, size = N, replace = FALSE)
    x =  cellLayerFromIndex(R, index)
    x$layer <- names(R)[x$layer]
    v =  layers_extractPoints(R, x)
    dplyr::bind_cols(x, v = v)
}


#' Select N random points from a faux multilayer Raster* object.
#'
#' Similar to \code{layers_randomPoints()} but here the user provides a vector of layer
#' identifiers (indicies, names, etc) that may or may not exist in the input raster, \code{R}.
#' Also, NA is returned as the value for each selected point.
#'
#' This function does not guarantee that returned locations do not coincide with
#' missing points.  The first layer in R is used with \code{dismo::randomPoints()}
#' as if the first layer is a mask so select M*N non-NA points.
#'
#' @export
#' @param R a Raster* object - treated as a single layers
#' @param layers a vector of layer IDs, by default 1,2, ..., nlayers(R) but can be anything
#' @param pts xyz values for presence points [lon, lat, layer] - these
#'  locations are avoided when sampling. Ignored if NULL.   This must have
#'  either set of the following columns.  Note that layer (or z)
#'  must be found in the \code{layer} argument.
#' \itemize{
#'      \item{'cell' and 'layer' (or 'z')}
#'      \item{'row', 'col' and 'layer' (or 'z')}
#'      \item{'lon', 'lat' and 'layer' (or 'z')}
#'      \item{'x', 'y' and 'layer' (or 'z')}
#'  }
#' @param N the number of points to select
#' @param M numeric, a multiplier to start with a selection of random points
#'          from which N are selected.
#' @return data.frame of sampled points with the following
#'  \itemize{
#'      \item{lon}
#'      \item{lat}
#'      \item{col}
#'      \item{row}
#'      \item{cell 2D matrix index}
#'      \item{index 3d array index}
#'      \item{layer}
#'      \item{value - everywhere NA to match \code{layers_randomPoints()}}
#'  }
faux_randomPoints <- function(R,
    layers = seq_len(raster::nlayers(R)),
    pts = NULL,
    N = 1000,
    M = 5){

    s           = raster_shape(R)
    s[['nl']]   = length(layers)
    s[['N']]    = s[['n']] * s[['nl']]

    xy  = dismo::randomPoints(R, n = N*M)
    z   = sample(seq_len(s[['nl']]), size = N*M, replace = TRUE)
    xyz = dplyr::tibble(x = xy[,'x'], y = xy[,'y'], layer = z)
    index = indexFromPts(R, xyz)
    if (!is.null(pts)){
        pindex = indexFromPts(R, pts)
        pindex = pindex[!is.na(pindex)]
        index  = index[!(index %in% pindex)]
    }
    index = sample(index, size = N, replace = FALSE)
    x =  cellLayerFromIndex(R, index, shape = s)
    x %>%
        dplyr::mutate(layer = layers[layer], value = NA)
}
