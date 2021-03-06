% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grd.R
\name{read_grd}
\alias{read_grd}
\title{Read the header portion of a raster file}
\usage{
read_grd(
  filename = file.path("/mnt/ecocast/coredata/obpg/gom/DAY/chlor_a",
    "A2014365.L3m_DAY_CHL_chlor_a_4km.grd"),
  form = c("list", "tibble")[2]
)
}
\arguments{
\item{filename}{on or more filenames}

\item{form}{return form can be 'list' or 'tibble' (default)}
}
\value{
one or more grd lists or the same as a tibble
\itemize{
     \item{dim vector of nrows, ncols, ncells, nlayers}
     \item{extent vector of xmin, xmax, ymin, ymax}
     \item{projection  string}
     \item{filename string}
 }
}
\description{
Read the header portion of a raster file
}
