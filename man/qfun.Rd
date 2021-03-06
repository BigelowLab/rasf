% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/queue.R
\name{qfun}
\alias{qfun}
\title{Transform a queue into a raster stack and compute the specified function
using \code{\link[raster]{calc}}}
\usage{
qfun(Q, fun = mean, ...)
}
\arguments{
\item{Q}{rdeque class object}

\item{fun}{the function to apply using \code{\link[raster]{calc}}}

\item{...}{further arguments for raster::calc where the function is mean}
}
\value{
RasterLayer class
}
\description{
Transform a queue into a raster stack and compute the specified function
using \code{\link[raster]{calc}}
}
