% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bb.R
\name{bb_to_bbox}
\alias{bb_to_bbox}
\title{Convert a 4-element bbox vector to a sf bbox object}
\usage{
bb_to_bbox(
  bb = c(-72, -63, 39, 46),
  crs = ifelse(use_wkt("sf"), "epsg:4326", "+init=epsg:4326")
)
}
\arguments{
\item{bb}{a 4-element numeric vector of left, right, bottom and top coordinates}

\item{crs}{character, the coordinate reference system}
}
\value{
sf bbox object
}
\description{
Convert a 4-element bbox vector to a sf bbox object
}
