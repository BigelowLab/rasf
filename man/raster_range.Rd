% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/raster.R
\name{raster_range}
\alias{raster_range}
\title{Computer the range of values for raster objects}
\usage{
raster_range(x, na.rm = TRUE, collapse = TRUE, ...)
}
\arguments{
\item{x}{a Raster* or SpatRaster class}

\item{na.rm}{logical, if TRUE remove NAs. For objects inheriting \code{BasicRaster}.
Ignored for class SpatRaster.}

\item{collapse}{logical, if FALSE then a matrix is returned with two columns
which are vectors of \code{min} and \code{max} values for each layer.
If TRUE then a a two element vector of  \code{min} and \code{max} for all layers
is returned}

\item{...}{other arguments for \code{\link[terra]{minmax}}}
}
\value{
either a matrix of min and max values or a two element vector of
   min and max
}
\description{
Computer the range of values for raster objects
}
