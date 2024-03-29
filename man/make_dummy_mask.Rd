% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example_data.R
\name{make_dummy_mask}
\alias{make_dummy_mask}
\title{Create a dummy raster mask for testing \code{make_ratser_lut}}
\usage{
make_dummy_mask(nc = 10, nr = 10, klass = c("RasterLayer", "SpatRaster")[2])
}
\arguments{
\item{nc}{integer, number of columns}

\item{nr}{integer, number of rows}

\item{klass}{character, one of \code{RasterLayer} of \code{SpatRaster}}
}
\value{
raster mask with NA values assigned to area to be masked
}
\description{
Create a dummy raster mask for testing \code{make_ratser_lut}
}
