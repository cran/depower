#' Simulate data from the null distribution
#'
#' A method to simulate data from the null distribution of
#' [depower::sim_log_lognormal()], [depower::sim_nb()], and [depower::sim_bnb()].
#'
#' @param data (depower)\cr
#'        The simulated data returned by [depower::sim_log_lognormal()],
#'        [depower::sim_nb()], or [depower::sim_bnb()].
#' @param .funs (named list of calls)\cr
#'        A list of named calls (functions) to be used for the tests in the
#'        power analysis.
#' @param ... Additional arguments passed to the function.
#'
#' @return A data frame
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # sim_null() examples
#' #----------------------------------------------------------------------------
#' library(depower)
#'
#' @noRd
sim_null <- function(
  data,
  .funs,
  ...
) {
  UseMethod("sim_null")
}

#' @export
sim_null.nb <- function(
  data,
  .funs,
  ...
) {
  data_and_fun_args(data, .funs) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      data_null = if (isTRUE(.data$distribution_test_stat == "simulated")) {
        sim_nb(
          n1 = .data$n1,
          n2 = .data$n2,
          mean1 = .data$mean1,
          ratio = if (is.na(.data$ratio_null)) 1 else .data$ratio_null,
          dispersion1 = .data$dispersion1,
          dispersion2 = if (isTRUE(.data$equal_dispersion)) {
            .data$dispersion1
          } else {
            .data$dispersion2
          },
          nsims = .data$nsims_test_stat
        ) |>
          getElement("data")
      } else {
        list(NULL)
      }
    ) |>
    dplyr::ungroup()
}

#' @export
sim_null.bnb <- function(
  data,
  .funs,
  ...
) {
  data_and_fun_args(data, .funs) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      data_null = if (isTRUE(.data$distribution_test_stat == "simulated")) {
        sim_bnb(
          n = .data$n1,
          mean1 = .data$mean1,
          ratio = if (is.na(.data$ratio_null)) 1 else .data$ratio_null,
          dispersion = .data$dispersion1,
          nsims = .data$nsims_test_stat
        ) |>
          getElement("data")
      } else {
        list(NULL)
      }
    ) |>
    dplyr::ungroup()
}

#' @export
sim_null.log_lognormal_one_sample <- function(
  data,
  .funs,
  ...
) {
  data_and_fun_args(data, .funs) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      data_null = if (isTRUE(.data$distribution_test_stat == "simulated")) {
        sim_log_lognormal(
          n1 = .data$n1,
          ratio = if (is.na(.data$mean_null)) 1 else exp(.data$mean_null),
          cv1 = .data$cv1,
          nsims = .data$nsims_test_stat
        ) |>
          getElement("data")
      } else {
        list(NULL)
      }
    ) |>
    dplyr::ungroup()
}

#' @export
sim_null.log_lognormal_independent_two_sample <- function(
  data,
  .funs,
  ...
) {
  data_and_fun_args(data, .funs) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      data_null = if (isTRUE(.data$distribution_test_stat == "simulated")) {
        sim_log_lognormal(
          n1 = .data$n1,
          n2 = .data$n2,
          ratio = if (is.na(.data$mean_null)) 1 else exp(.data$mean_null),
          cv1 = .data$cv1,
          cv2 = .data$cv2,
          cor = .data$cor,
          nsims = .data$nsims_test_stat
        ) |>
          getElement("data")
      } else {
        list(NULL)
      }
    ) |>
    dplyr::ungroup()
}

#' @export
sim_null.log_lognormal_dependent_two_sample <- function(
  data,
  .funs,
  ...
) {
  data_and_fun_args(data, .funs) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      data_null = if (isTRUE(.data$distribution_test_stat == "simulated")) {
        sim_log_lognormal(
          n1 = .data$n1,
          n2 = .data$n2,
          ratio = if (is.na(.data$mean_null)) 1 else exp(.data$mean_null),
          cv1 = .data$cv1,
          cv2 = .data$cv2,
          cor = .data$cor,
          nsims = .data$nsims_test_stat
        ) |>
          getElement("data")
      } else {
        list(NULL)
      }
    ) |>
    dplyr::ungroup()
}

#' @export
sim_null.log_lognormal_mixed_two_sample <- function(
  data,
  .funs,
  ...
) {
  data_and_fun_args(data, .funs) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      data_null = if (isTRUE(.data$distribution_test_stat == "simulated")) {
        sim_log_lognormal(
          n1 = .data$n1,
          n2 = .data$n2,
          ratio = if (is.na(.data$mean_null)) 1 else exp(.data$mean_null),
          cv1 = .data$cv1,
          cv2 = .data$cv2,
          cor = .data$cor,
          nsims = .data$nsims_test_stat
        ) |>
          getElement("data")
      } else {
        list(NULL)
      }
    ) |>
    dplyr::ungroup()
}

#-------------------------------------------------------------------------------
# Helpers
#-------------------------------------------------------------------------------
fun_args <- function(.funs) {
  lapply(
    X = .funs,
    FUN = function(x) {
      dist <- eval(x$distribution)
      list(
        distribution = dist$distribution %||% NA,
        nsims = dist$nsims %||% NA,
        ncores = dist$ncores %||% NA,
        alternative = x$alternative %||% NA,
        equal_dispersion = x$equal_dispersion %||% NA,
        ratio_null = x$ratio_null %||% NA,
        mean_null = x$mean_null %||% NA
      )
    }
  )
}

data_and_fun_args <- function(data, .funs) {
  .fun_args <- fun_args(.funs)

  data |>
    dplyr::select(-dplyr::any_of("data")) |>
    dplyr::cross_join(
      y = list2DF(
        list(
          test = setNames(names(.funs), as.character(.funs)),
          distribution_test_stat = unlist(
            lapply(.fun_args, function(x) x$distribution),
            use.names = FALSE
          ),
          nsims_test_stat = unlist(
            lapply(.fun_args, function(x) x$nsims),
            use.names = FALSE
          ),
          ncores_test_stat = unlist(
            lapply(.fun_args, function(x) x$ncores),
            use.names = FALSE
          ),
          alternative = unlist(
            lapply(.fun_args, function(x) x$alternative),
            use.names = FALSE
          ),
          equal_dispersion = unlist(
            lapply(.fun_args, function(x) x$equal_dispersion),
            use.names = FALSE
          ),
          ratio_null = unlist(
            lapply(.fun_args, function(x) x$ratio_null),
            use.names = FALSE
          ),
          mean_null = unlist(
            lapply(.fun_args, function(x) x$mean_null),
            use.names = FALSE
          )
        )
      )
    )
}
