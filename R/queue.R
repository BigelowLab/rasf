#' Add an element to the front of a queue
#'
#' @seealso \code{\link[rstackdeque]{insert_front}}
#' @export
#' @param Q rdeque class object
#' @param x an element to add to the front of the queue
#' @return the updated queue
qpush    <- function(Q, x) {
  rstackdeque::insert_front(Q, x)
}

#' Add an element to the front of a queue
#'
#' @seealso \code{\link[rstackdeque]{without_back}}
#' @export
#' @param Q rdeque class object
#' @return the removed queue element
qpop <- function(Q){
  invisible(rstackdeque::without_back(Q))
}

#' Add an element to the front of a queue
#'
#' @export
#' @param Q rdeque class object
#' @param x an element to add to the queue, the forst on is dropped
#' @return the updated queue
qadvance <- function(Q, x){
  qpush(qpop(Q), x)
}

#' Transform a queue into a raster stack
#'
#' @export
#' @param Q rdeque class object
#' @return RasterStack class
qstack   <- function(Q) {
  raster::stack(Q)
}

#' Transform a queue into a raster stack and compute the mean
#'
#' @export
#' @param Q rdeque class object
#' @param ... furthere arguments for \code{\link[raster]{calc}} where the
#'        function is mean
#' @return RasterLayer class
qavg     <- function(Q, ...) {
  qfun(Q, fun = mean, ...)
}

#' Transform a queue into a raster stack and compute the specified function
#' using \code{\link[raster]{calc}}
#'
#' @export
#' @param Q rdeque class object
#' @param fun the function to apply using \code{\link[raster]{calc}}
#' @param ... further arguments for raster::calc where the function is mean
#' @return RasterLayer class
qfun     <- function(Q, fun = mean, ...) {
  raster::calc(raster::stack(as.list(Q)), fun, ...)
}

#' Create a queue
#'
#' @export
#' @param x a raster to apply as the first element
#' @return rdeque class
qinit <- function(x = NULL){
  Q <- rstackdeque::rdeque()
  if (!is.null(x)) Q <- qpush(Q, x)
  Q
}
