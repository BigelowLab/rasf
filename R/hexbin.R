#' Retrieve of set of hex-binned polygons for the given xyz data
#'
#' @export
#' @seealso \href{https://www.rdocumentation.org/packages/hexbin/versions/1.29.0/topics/hexbin}{hexbin docs}
#' @seealso \href{https://hexnet.org/content/hexagonal-geometry}{nice overview}
#' @param x a vector of lon-locations
#' @param y a vector of lat-locations
#' @param z a vector of values associated with each location
#' @param fun character the function used to process the cell values
#'    currently only 'mean', 'median', or 'count' (default)
#' @param crs the projection string for final data
#' @param utm_crs intermediary UTM projection string
#' @param na.rm logical passed through to user's selected function
#' @param ... further arguments for \code{hexbin::hexbin()}
#' @return a POLYGON sf data frame
st_hexbin_utm <- function(x, y, z,
                    fun = c('mean', 'count', 'median')[2],
                    crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
                    utm_crs = "+proj=utm +zone=19 ellps=WGS84",
                    na.rm = TRUE,
                    ...){



    pts <- sf::st_as_sf(x= dplyr::tibble(lon = x, lat = y),  coords = c("lon", "lat"), crs = crs)
    utm <- sf::st_transform(pts, crs = utm_crs)
    xy <- sf::st_coordinates(utm)
    hb  <- hexbin::hexbin(xy[,1], xy[,2], IDs = TRUE, ...)
    # compute the value
    if (fun %in% c('mean', 'median')){
        cID     <- hb@cID
        cell    <- hb@cell
        idx     <- split(seq_along(cID), cID)
        v       <- rep(0, hb@ncells)
        for (i in seq_along(cell)) {
            ix <- idx[[as.character(cell[i])]]
            v[i] <- switch(fun,
                           "mean" = mean(z[ix], na.rm = na.rm),
                           "median" = median(z[ix], na.rm = na.rm)
                            )
        }
    } else {
        v <- hb@count
    }

    rx <- diff(hb@xbnds)/hb@xbins/2
    ry <- diff(hb@ybnds)/hb@xbins/2
    R <- 2/sqrt(3)*rx
    cos30 <- cos(30*pi/180)
    sin30 <- sin(30*pi/180)

    px <- R * c(0, cos30,  cos30,  0, -cos30, -cos30, 0)
    py <- R * c(1, sin30, -sin30, -1, -sin30,  sin30, 1)

    xy <- hexbin::hcell2xy(hb)
    xc <- xy$x
    yc <- xy$y
    cell <- hb@cell

    pp <- lapply(seq_along(xc),
        function(i){
            sf::st_polygon(list(cbind(xc[i] + px, yc[i] + py)))
        })
    sf::st_transform(sf::st_sf(value = v, cell = pp, crs = utm_crs), crs = crs)
}




#' Retrieve of set of hex-binned polygons for the given xyz data
#'
#' @export
#' @seealso \href{https://www.rdocumentation.org/packages/hexbin/versions/1.29.0/topics/hexbin}{hexbin docs}
#' @seealso \href{https://hexnet.org/content/hexagonal-geometry}{nice overview}
#' @param x a vector of x-locations
#' @param y a vector of y-locations
#' @param z a vector of values associated with each location
#' @param fun character the function used to process the cell values
#'    currently only 'mean', 'median', or 'count' (default)
#' @param crs the projection string
#' @param na.rm logical passed through to user's selected function
#' @param ... further arguments for \code{hexbin::hexbin()}
#' @return a POLYGON ssf data frame
st_hexbin <- function(x, y, z,
                    fun = c('mean', 'count', 'median')[2],
                    crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
                    na.rm = TRUE,
                    ...){

    if (missing(x)){
        xyz <- readr::read_csv("/mnt/ecocast/projectdata/moosecrash/versions/v0/v0.000/xcast/278/2018-10-05-06.csv.gz")
        x <- xyz$lon
        y <- xyz$lat
        z <- xyz$xcast
        fun = 'mean'
        crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
        na.rm = TRUE
    }

    # compute the width in meters of one degree of longitude which depends upon latitude
    # https://gis.stackexchange.com/questions/142326/calculating-longitude-length-in-miles
    # d_lon <- cos(lat) * del_circ
    lon_width_at_lat <- function(lat, del_circ = 40075017/360){
        cos(lat * pi/180) * del_circ
    }
    # https://en.wikipedia.org/wiki/Longitude#Length_of_a_degree_of_longitude
    lon_length_at_lat <- function(lat, a = 6378137.0, b = 6356752.3142){
        e = (a^2 - b^2)/a^2
        theta = lat * pi/180
        (pi*a*cos(theta))/(180*sqrt(1-e^2*sin(theta)^2))
    }
    # https://en.wikipedia.org/wiki/Latitude#Length_of_a_degree_of_latitude
    lat_length_at_lat <- function(lat, p = c(a=111132.954, b = 559.822, c = 1.175  ) ){
        theta = lat * pi/180
        p[['a']] - p[['b']]*cos(2*theta) + p[['c']]*cos(4*theta)
    }

    stopifnot(all.equal(length(x), length(y), length(z)))

    fun <- tolower(fun[1])

    # make the hexbin object
    hb      <- hexbin::hexbin(x, y, IDs = TRUE, ...)

    # compute the value
    if (fun %in% c('mean', 'median')){
        cID     <- hb@cID
        cell    <- hb@cell
        idx     <- split(seq_along(cID), cID)
        v       <- rep(0, hb@ncells)
        for (i in seq_along(cell)) {
            ix <- idx[[as.character(cell[i])]]
            v[i] <- switch(fun,
                           "mean" = mean(z[ix], na.rm = na.rm),
                           "median" = median(z[ix], na.rm = na.rm)
                            )
        }
    } else {
        v <- hb@count
    }


    f <- 2/sqrt(3)
    theta <- 360/6
    bearings <- c(theta * 0:5, 0) # N, NE, SE, S, SW, NW, N

    # now set up the polygons
    # get the cell centers (not centers of mass as in xcm, ycm)
    xy <- hexbin::hcell2xy(hb)
    xc <- xy$x
    yc <- xy$y
    cell <- hb@cell
    # from hexbin::hexbin() docs
    # hexagon inner diameter = diff(xbnds)/xbins in x units
    # in this case degrees which we convert to meters depending upon latitude
    w <- diff(hb@xbnds)/hb@xbins/2
    dx <- lon_length_at_lat(yc)
    dy <- lat_length_at_lat(yc)
    r <- w * dx
    RNS  <- w * dy
    RB <- r*f
    R <- cbind(RNS, RB, RB, RNS, RB, RB, RNS)  # N, NE, SE, S, SW, NW, N



    # destPoint(p, b, d)
    pp <- lapply(seq_along(xc),
        function(i) {
            p <- sapply(seq_along(bearings),
                function(j){
                    geosphere::destPoint(c(xc[i], yc[i]), bearings[j], R[i,j])
            })
            sf::st_polygon(list(t(p)))
        })

    sf::st_sf(value = v, cell = pp, crs = crs)

}
