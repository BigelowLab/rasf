# #' Perform binary dilation on a binary raster
# #'
# #' @export
# #' @param x binary raster where the dilated portion has value of 1, all other
# #'    cells have value 0
# #' @param width numeric odd numbered kernal width
# #' @param height numeric odd numbered kernal width
# #' @param pad logical, if TRUE pad the input to resolve edge effects
# #' @param padValue numeric value of pad
# #' @return raster layer
# raster_dilate <- function(x, width = 3, height = width, pad = TRUE, padvalue = 0){
#
#     k = matrix(1/(width[1] * height[1]), ncol = width[1], nrow = height[1])
#
#     ones        = x[] == 1
#     zeroes      = !ones
#
#     x2          = raster::focal(x, f, pad = pad, padValue = padValue)
#
#     x2[ones]    = 1
#     swath       = x2[]
#     swath       = swath > 0 & swath < 1
#
#     x3          = x
#
#
#
# library(raster)
# nr = 30
# nc = 30
# m = matrix(seq_len(nc*nr), ncol = nc, nrow = nr)
# r = raster(upper.tri(m, diag = TRUE))
# ocean   = r[] == 0
# land    = !ocean
#
# fr = 5
# fc = 5
# f = matrix(1/(fr*fc), nrow = fr, ncol = fr)
#
# r2 = focal(r,f, pad = TRUE, padValue = 0)
#
# r2[land]    = 1
# swath       = r2[]
# swath       = swath > 0 & swath < 1
# r3          = r
# r3[]        = 1
# r3[swath]   = 0
#
# s           = stack(r, r2, r3)
# names(s)    = c("mask", "focal", "swath")
# plot(s)
#
#
#
# }



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
read_grd <- function(filename = '/mnt/ecocast/coredata/obpg/gom/DAY/chlor_a/A2014365.L3m_DAY_CHL_chlor_a_4km.grd',
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

#' Convert a index (1... ncells*nlayers) into a multilayer raster to cell and
#'    layer coordinates
#'
#' @export
#' @param index vector of multilayer index coordinates
#' @param R brick or stack
#' @return a tibble of index, cell, col, row, x, y, and layer
cellLayerFromIndex <- function(R, index){

    nr = raster::nrow(R)
    nc = raster::ncol(R)
    n  = raster::ncell(R)
    nl = raster::nlayers(R)
    N  = nc*nr*nl
    layer  = ((index-1) %/% n) + 1
    cell = index - ((layer - 1) * n)
    ix = ((cell-1) %% nc)  + 1
    iy = floor((cell - 1) / nc) + 1
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
    nc <- as.integer(raster::ncell(R))
    nl <- as.integer(raster::nlayers(R))
    ny <- raster::nrow(R)
    nx <- raster::ncol(R)

    pts <- vet_pts(pts)
    # rows
    x   <- raster::colFromX(R, pts$x)
    # cols
    y   <- raster::rowFromY(R, pts$y)
    z   <- pts$layer
    if (is.character(layer)) layer <- match(layer, names(R))
    # row + number of rows * rowsize + number of layers * layersize
    index <- x + (y-1)*nx  + (z-1)*nc
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
    nc <- as.integer(raster::ncell(R))
    nl <- as.integer(raster::nlayers(R))
    ny <- raster::nrow(R)
    nx <- raster::ncol(R)

    if (!(is.data.frame(pts) || is.matrix(pts)))
        stop("pts must be data.frame or matrix")
    if (inherits(pts, "matrix") || inherits(pts, 'tbl_df'))
        pts <- as.data.frame(pts, stringsAsFactors = FALSE)

    nm <- names(pts)
    iz <- which(nm %in% c("layer", "z"))[1]
    if (length(iz) == 0) stop("pts must have 'layer' or 'z' column")
    layer <- pts[,iz]
    if (is.character(layer)) layer <- match(layer, names(R))

    ix <- which(nm %in% c("x", "lon"))[1]
    if (length(ix) == 0) stop("pts must have 'x' or 'lon' column")

    iy <- which(nm %in% c("y", "lat"))[1]
    if (length(iy) == 0) stop("pts must have 'y' or 'lat' column")

    flayer  <- factor(layer)
    pp      <- split(pts, flayer)
    ii      <- as.numeric(levels(flayer))
    vv      <- sapply(ii,
        function(i){
            raster::extract(R[[i]], pp[[i]][,c(ix, iy)])
        }
    )
    v       <- unsplit(vv, flayer, drop = FALSE)

    #cell <- cellFromPts(R, pts)

    #index <- cell + (layer-1) * nc

    #raster::getValues(R)[index]
    unname(v)
}

#' Select N non-NA random points from a mulitlayer Raster* object
#'
#' @export
#' @param R a multilayer Raster* object
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
#' @param M numeric, a multiplier to start with a selection of N*M random points
#' from which N non-NA and non-pts are selected
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
layers_randomPoints <- function(R, pts = NULL, N = 1000, M = 5){

    if (!inherits(R, 'BasicRaster')) stop("Input R must be a Raster* class")
    nc <- raster::ncell(R)
    nl <- raster::nlayers(R)
    ny <- raster::nrow(R)
    nx <- raster::ncol(R)

    if (N > (nc * nl))
        stop("number of random points requested, N, exceeds number of available cells")
    if ((N*M) > (nc * nl))
        stop("number of pooled random points requested, N*M, exceeds available number of cells")

    index   = sample(nc * nl, N*M)
    xyz     = cellLayerFromIndex(R, index)
    v       = layers_extractPoints(R, xyz)
    vna     = is.na(v)

    if (!is.null(pts)){
        pts = vet_pts(pts)
        pts_idx = indexFromPts(R, pts)
        ix <- match(pts_idx, index)
        ixna <- is.na(ix)
        if (!all(ixna)){
            index[ix[!ixna]] <- NA
        }
    }
    if (any(vna)) index[vna] <- NA
    index <- index[!is.na(index)]
    if (length(index) < N)
        stop("N exceeds number of available values to select")
    index   = sample(index, N)
    xyz     = cellLayerFromIndex(R, index)
    v       = layers_extractPoints(R, xyz)
    x       = dplyr::tibble(
                lon    = xyz$x,
                lat    = xyz$y,
                row    = xyz$row,
                col    = xyz$col,
                cell   = xyz$cell,
                index  = index,
                layer  = xyz$layer,
                value  = v)
    invisible(x)

}

# layers_randomPoints_deprecated <- function(R, pts = NULL, N = 1000, M = 5){
#
#     if (!inherits(R, 'BasicRaster')) stop("Input R must be a Raster* class")
#     if (!is.null(pts)){
#         if (!(is.data.frame(pts) || is.matrix(pts)))
#             stop("pts must be data.frame or matrix")
#         if (inherits(pts, 'tbl_df')) pts <- as.data.frame(pts, stringsAsFactors = FALSE)
#         if (ncol(pts) < 3) stop("pts must have at least three columns")
#         if (!any(c("lon", "x", "col") %in% names(pts)))
#             stop("pts must have 'col', 'lon' or 'x' column")
#         if (!any(c("lat", "y", "row") %in% names(pts)))
#             stop("pts must have 'row', 'lat' or 'y' column")
#         iz <- which(names(pts) %in% c("layer", "z"))[1]
#         if (length(iz) == 0) stop("pts must have 'layer' or 'z' column")
#     }
#     # get all of the data as a vector - storage is 1,2,3 across rows
#     # starting form upper left
#     v <- raster::getValues(R)
#     nc <- as.integer(raster::ncell(R))
#     nl <- as.integer(raster::nlayers(R))
#     ny <- raster::nrow(R)
#     nx <- raster::ncol(R)
#
#
#     # if pts are present, then we flag those in vector 'v'
#     if (!is.null(pts)){
#         if (is.character(pts[,iz])){
#             z <- match(pts[,iz], names(R))
#         } else {
#             z <- match(pts[,iz], 1:nl)
#         }
#         #col <- colFromX(R, pts[, ix])
#         #row <- rowFromY(R, pts[, iy])
#         #cell <- cellFromRowCol(R, row, col)
#         cell <- cellFromPts(R, pts)
#         # this is the number of layer steps 'extra'
#         bump <- 0:(nl-1L) * nc
#         index <- cell + bump[z]
#         v[index] <- NA
#     }
#     # create a dummy index
#     aix <- seq.int(from = 1L, to = (nc * nl))
#     # determine where NAs occur
#     nna <- !is.na(v)
#     if (sum(nna) < N) stop("N exceeds number of available values to select")
#     # sample the *indices*
#     s <- sort(sample(aix[nna], N))
#     # convert to cells, layers, rows, columns, and xy
#     cell <- ((s-1L) %% nc) + 1L
#     layer <- ((s-1L) %/% nc) + 1L
#     rc <- raster::rowColFromCell(R, cell)
#     xy <- raster::xyFromCell(R, cell)
#
#     x <- data.frame(
#         lon = raster::xFromCol(R, rc[,'col']),
#         lat = raster::yFromRow(R, rc[,'row']),
#         row = rc[,'row'],
#         col = rc[,'col'],
#         cell = cell,
#         layer = layer,
#         value = v[cell + (layer-1)*nc],
#         stringsAsFactors = TRUE)
#     invisible(x)
# }
