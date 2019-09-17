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
#' @return a matrix of length(x) rows with columns nrow, ncol, ncell, nlayer
grd_dim <- function(x){
  x <- t(sapply(x, '[[', 'dim'))
  colnames(x) = c("nrow", "ncol", "ncell", "nlayer")
  x
}
