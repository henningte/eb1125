#' Saving and loading tibbles with rvars columns
#'
#' See https://github.com/stan-dev/posterior/issues/307
#'
#' @param object Object to save.
#'
#' @param file Character value. Path to the file.
#'
#' @export
saveRDS_rvars <- function(object, file) {

  saveRDS(
    object = object,
    file = file,
    refhook = \(x) if (any(c("vec_proxy", "vec_proxy_equal") %in% names(x))) ""
  )

}

#' @export
readRDS_rvars <- function(file) {

  readRDS(
    file = file,
    refhook = \(x) new.env()
  )

}
