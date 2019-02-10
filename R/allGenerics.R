#' @export
setGeneric("biosign", function(x, ...) standardGeneric("biosign"))

#' @rdname getAccuracyMN
#' @export
setGeneric("getAccuracyMN", function(object, ...) {standardGeneric("getAccuracyMN")})

#' @rdname getSignatureLs
#' @export
setGeneric("getSignatureLs", function(object, tierC = c("S", "AS")[1], ...) {standardGeneric("getSignatureLs")})
