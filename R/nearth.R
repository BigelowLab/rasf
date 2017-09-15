#' Various bounding boxes for analyses or forecasts
#'
#' Note - analyses are LCC coordinates while forecasts are longlat 
#'
#' @export
#' @param where character 'native' (aka 'all'), 'nwa', 'neac' etc
#' @return 4 element vector of [left, right, bottom, top]
celmap_bb <- function(
    where = c('all', 'maine', 'gom', 'nwa', 'neac', 'native')[2], 
    pad = list(a = c(0,0), b = c(0.1, 0.1))[[1]],
    form = c("numeric", "extent", "sp")[1]){

        
    bb = switch(tolower(where[1]),
        'maine' =   c(-71.1, -67, 43, 47.5),
        'gom'   =   c(-72,-63,39,46),
        'nwa'   =   c(-77.1199677903564, -51.6967983980479, 
                      37.9115974449456, 56.6768926463509),
        'neac'  =   c(-74, -59.75, 41, 48.15),
                    c(-180, 180, -90, 90))
    if (!is.null(pad) && !is.na(pad)){
        if (length(pad) == 0) pad = c(pad, pad)
        bb = bb + c(-pad[1], pad[1], -pad[2], pad[2])
    }
    
    switch(tolower(form[1]),
        'extent' = raster::extent(bb),
        'sp' = nearth::bbox_to_SpatialPolygons(bb),
        bb)
}


#' Get a Spatial vector object by region name and vector set.
#' 
#' Currently draws from Natural Earth
#'
#' @export
#' @param where the name of the region
#' @param ... further arguments for \code{celmap_bb()}
#' @return a Spatial* object
get_vectors = function(where = c("World", "Maine", "GOM", "NWA", "NEAC")[1], 
    what = c('coast50m', "coast10m", 'boundary50m', 'boundary10m')[1], ...){
    
    what = switch(tolower(what[1]),
        'coast50m' = 'ne_50m_coastline',
        'boundary50m' = 'ne_50m_admin_1_states_provinces_shp',
        'coast10m' = 'ne_10m_coastline',
        'boundary10m' = 'ne_10m_admin_1_states_provinces_shp',
        tolower(what[1]))
        
    bb = celmap_bb(where[1], ...)
    
    nearth::read_nearth(what,what = 'vector', bb = bb)
}