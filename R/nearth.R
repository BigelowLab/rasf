#' Get a Spatial vector object by state or province name
#'
#' 'Quebec' can be used in lieu of "Qu\xe9bec"
#' @export
#' @param where names the states and or provinces, or group names such as 'GOM'
#' @param what the dataset to quuery (boundary50m - the default - or boundary10m)
#' @return a Spatial* object
get_boundaries = function(
    where = c("Maine", "New Brunswick", "Quebec", "New Hampshire"), 
    what = c('boundary50m', 'boundary10m')[1]){
    
    
    wh = switch(tolower(where[1]),
        'gom' = c('Rhode Island', 'Massachusetts', 'New Hampshire',
            'Vermont', 'Quebec', 'New Brunswick', 'Nova Scotia' , 'Maine'),
        where)
        
    ix = tolower(wh) %in% 'quebec'
    if (any(ix)) wh[ix] = "Qu\xe9bec"
    
    
    what = switch(tolower(what[1]),
        'boundary50m' = 'ne_50m_admin_1_states_provinces_shp',
        'boundary10m' = 'ne_10m_admin_1_states_provinces_shp',
        tolower(what[1]))
    
    x = nearth::read_nearth(what, what = 'vector')[[1]]
    
    x[x@data$name %in% where,]
}
     
#' Get a Spatial vector object by region name and vector set.
#' 
#' Currently draws from Natural Earth.  Note that this clips the region
#' so edge effects appear - especially for boundaries.  See \code{get_boundary()}
#'
#' @export
#' @param where the name of the region
#' @param what the type of boundary (10m, 50m etc)
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
        
    bb = get_bb(where[1], ...)
    
    nearth::read_nearth(what,what = 'vector', bb = bb)[[1]]
}