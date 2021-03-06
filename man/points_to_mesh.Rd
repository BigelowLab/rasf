% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sf.R
\name{points_to_mesh}
\alias{points_to_mesh}
\title{Convert sfc_POINTS to a mesh of polygons (delaunay triangulation)}
\usage{
points_to_mesh(x, varname = NULL, fun = mean, ...)
}
\arguments{
\item{x}{sf object of type sfc_POINT}

\item{varname}{character one or more variable names to transfer to the polygon}

\item{fun}{function to compute value the variables specifed by varname}

\item{...}{other arguments for \code{fun}}
}
\value{
sfc_POLYGON table
}
\description{
Convert sfc_POINTS to a mesh of polygons (delaunay triangulation)
}
