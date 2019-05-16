#' Various bounding boxes for analyses or forecasts
#'
#'
#' @export
#' @param where character 'native' (aka 'all' or 'world'), 'nwa', 'neac' etc
#' \itemize{
#'   \item{'maine' \code{-71.1, -67, 43, 47.5}}
#'   \item{'gom' Gulf of Maine \code{-72, -63, 39, 46}}
#'   \item{'nwa' Northwest Atlantic \code{-77.0, -51.5, 37.9, 56.7}}
#'   \item{'neac' New England Atlantic Canada \code{-74, -59.75, 41, 48.15}}
#'   \item{'liac' Long Island Atlantic Canada \code{-74, -59.75, 37.9, 48.15}}
#'   \item{'gosl' Gulf of St Lawrence \code{-67,-56.5,44.4,50.5}}
#'   \item{any other \code{-180, 180, -90, 90}}
#' }
#' @param pad numeric one or two element pad to add to a bounding box,
#'     default is NULL to skip
#' @param form character, either "numeric", "extent", "sf", "sp",
#'    "top-left-bottom-right" or "left-top-right-bottom"
#' @param ... further arguments - likely proj_string for sp and sf outputs
#' @return 4 element vector of [left, right, bottom, top], extent or SpatialPolygons object
get_bb <- function(
    where = c('all', 'maine', 'gom', 'nwa', 'neac', 'native')[2],
    pad = NULL,
    form = c("numeric", "extent", "sp", "sf", "left-top-right-bottom",
             "top-left-bottom-right")[1],
    ...){


    bb = switch(tolower(where[1]),
        'maine' =   c(-71.1, -67, 43, 47.5),
        'gom'   =   c(-72,-63,39,46),
        #'nwa'   =   c(-77.1199677903564, -51.6967983980479,
        #              37.9115974449456, 56.6768926463509),
        'nwa'   =   c(-77.0, -51.5, 37.9, 56.7),
        'neac'  =   c(-74, -59.75, 41, 48.15),
        'liac' =    c(-74, -59.75, 37.9, 48.15),
        'gosl'  =   c(-67, -56.5, 44.4, 50.5),
                    c(-180, 180, -90, 90))

    if (!is.null(pad) && !is.na(pad)){
        if (length(pad) == 1) pad = c(pad, pad)
        bb = bb + c(-pad[1], pad[1], -pad[2], pad[2])
    }

    switch(tolower(form[1]),
        'extent'    = raster::extent(bb),
        'sp'        = bbox_to_SpatialPolygons(bb, ...),
        'sf'        = bbox_to_sf(bb, ...),
        "left-top-right-bottom" = bb[c(1,4,2,3)],
        "top-left-bottom-right" = bb[c(4,1,3,2)],
                      bb)
}
