#' @examples
#' #----------------------------------------------------------------------------
#' # asymptotic() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' set.seed(1234)
#' data <- sim_nb(
#'   n1 = 60,
#'   n2 = 40,
#'   mean1 = 10,
#'   ratio = 1.5,
#'   dispersion1 = 2,
#'   dispersion2 = 8
#' )
#'
#' data |>
#'   wald_test_nb(distribution = asymptotic())
#'
#' @name distribution
NULL

#' @export
#' @rdname distribution
asymptotic <- function() {
  list(
    distribution = "asymptotic",
    method = "asymptotic"
  )
}

check_asymptotic <- function(x) {
  if (!identical(x, asymptotic())) {
    stop(
      "Check argument 'distribution'. You must use 'distribution = depower::asymptotic()'."
    )
  }

  invisible(x)
}
