#' Retrieve of set of hex-binned polygons for the given xyz data
#'
#' @seealso \href{https://www.rdocumentation.org/packages/hexbin/versions/1.29.0/topics/hexbin}{hexbin docs}
#' @seealso \href{https://hexnet.org/content/hexagonal-geometry}{nice overview}
#' @param x a vector of x-locations
#' @param y a vector of y-locations
#' @param z a vector of values associated with each location
#' @param fun character the function used to process the cell values
#'    currently only 'mean', 'median', or 'count' (default)
#' @param crs the projection string
#' @param ... further arguments for \code{hexbin::hexbin()}
#' @return a sf data frame
st_hexbin <- function(x, y, z,
                    fun = c('mean', 'count', 'median')[2],
                    crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
                    ...){

    fun <- tolower(fun[1])

    # make the hexbin object
    hb      <- hexbin::hexbin(x, y, IDs = TRUE)

    # compute the value
    if (fun[1] %in% c('mean', 'median')){
        cID     <- hb@cID
        cell    <- hb@cell
        idx     <- split(seq_along(cID), cID)
        v       <- rep(0, hb@ncells)
        for (i in seq_along(cell)) {
            ix <- idx[[as.character(cell[i])]]
            v[i] <- switch(fun,
                           "mean" = mean(xyz$z[ix], na.rm = na.rm),
                           "median" = percentile(xyz$z[ix], probs = 0.5, na.rm = na.rm, names = FALSE)
                            )
        }
    } else {
        v <- hb@count
    }

    # now set up the polygons
    xc <- hb@xcm
    yc <- hb@ycm
    cell <- hb@cell
    # from hexbin::hexbin() docs
    # hexagon inner diameter =diff(xbnds)/xbins in x units
    #
    r <- diff(hb@xbnds)/hb@xbins/2
    r2 <- r/2
    R <- sqrt(3)/2*r
    px <- c(r,r,0,-r,-r,0,r)
    py <- c(r2,-r2,-R,-r2,r2,R,r2)
    pp <- lapply(seq_along(xc),
        function(i) sf::st_polygon(list(cbind(xc[i] + px, yc[i] + py)))
        )

    sf::st_sf(value = v, pp, crs = crs)
}
