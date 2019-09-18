#' Transform locations in a table from one projection to another
#'
#' Elements of the input that are NA are returned as NA
#'
#' @export
#' @param x table(data.frame or tibble) with columns specified by \code{from_names}
#' @param from_crs projection of the source coordinates - see \code{\link[sf]{st_crs}}
#' @param to_crs projection of the source coordinates - see \code{\link[sf]{st_crs}}
#' @param from_names string, the names of the x and y coords to project from, ['lon', 'lat'] default.
#' @param to_names string, the names of the x and y coords to project into, ['x', 'y'] default.
#'        If from_names is the same as to_names then replacement occurs, otherwise two
#'        new columns \code{to_names} are appended to the table.
#' @return updated table
#' @examples
#' \dontrun{
#' x <- data.frame(lon = c(-72, -63), lat = c(39, 46))
#' (y <- project_table(x, to_names = c("false_easting", "false_northing")))
#' # lon lat false_easting false_northing
#' # 1 -72  39      240199.8        4321059
#' # 2 -63  46      964564.1        5111577
#' }
project_table <- function(x = data.frame(lon = c(-72, -63), lat = c(39, 46)),
                           from_crs = "+init=epsg:4326",
                           to_crs = "+init=epsg:32619",
                           from_names = c('lon', 'lat'),
                           to_names = c("x", "y")){

  ll <- x %>% dplyr::select(from_names)
  ix <- apply(ll, 1, function(x) any(is.na(x)) )
  #input <- ll %>% dplyr::filter(!ix)
  output = sf::st_as_sf(ll %>% dplyr::filter(!ix),
                   coords = from_names,
                   crs = from_crs,
                   agr = "identity") %>%
          sf::st_transform(crs = to_crs) %>%
          sf::st_coordinates() %>%
          dplyr::as_tibble()
  #sp::coordinates(input) <- from_names
  #sp::proj4string(input) <- from_crs
  #output <- sp::coordinates(sp::spTransform(input, to_crs))
  ll[!ix,] <- output
  #ll[!ix,] <- tibble::as_tibble(output[,1])
  colnames(ll) <- to_names
  if (tibble::has_name(x, to_names[1])) {
    x[[to_names[[1]]]] <- ll[[to_names[1]]]
  } else {
    x <- x %>% tibble::add_column(!!to_names[1] := ll[[to_names[1]]])
  }
  if (tibble::has_name(x, to_names[2])) {
    x[[to_names[[2]]]] <- ll[[to_names[2]]]
  } else {
    x <- x %>% tibble::add_column(!!to_names[2] := ll[[to_names[2]]])
  }
  x
}
