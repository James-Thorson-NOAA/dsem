
#' @title
#' Family for data that are known without error
#'
#' @description
#' Allows using \code{family = fixed()} to specify data that have no measurement error
#'
#' @export
fixed <- function() {
  link1 = "identity"
  l1 <- substitute(link1)
  if (!is.character(l1)) l1 <- deparse(l1)
  .type <- "fixed"
  clean_name <- "fixed"
  structure(
    list(
      link = l1,
      type = .type,
      family = "fixed",
      clean_name = clean_name
    ),
    class = "family"
  )
}

#' @title
#' Gaussian with known standard deviation for measurement errors
#'
#' @description
#' Allows using \code{family = gaussian_fixed_sd(link, sd)} to specify data that have known measurement error
#'
#' @param link Link function
#' @param sd vector of known standard deviations (repeated for length of time-series, so be careful)
#'
#' @export
gaussian_fixed_sd <- function( link, sd ) {
  if( missing(link) ){
    link = "identity"
  }
  l1 <- substitute(link)
  if (!is.character(l1)) l1 <- deparse(l1)
  .type <- "gaussian_fixed_sd"
  clean_name <- "gaussian_fixed_sd"

  structure(
    list(
      link = l1,
      fixed_sd = sd,
      type = .type,
      family = "gaussian_fixed_sd",
      clean_name = clean_name
    ),
    class = "family"
  )
}

#' @title
#' Lognormal distribution
#'
#' @description
#' Lognormal distribution
#'
#' @details
#' Central tendancy meanlog is the median of log-response
#' i.e., does not include bias correction term $sigma^2/2$
#'
#' @param link Link function
#' @param sd vector of known standard deviations (repeated for length of time-series, so be careful)
#'
#' @export
lognormal <- function( link ) {
  if( missing(link) ){
    link = "log"
  }
  l1 <- substitute(link)
  if (!is.character(l1)) l1 <- deparse(l1)

  structure(
    list(
      link = l1,
      fixed_sd = sd,
      type = "lognormal",
      family = "lognormal",
      clean_name = "lognormal"
    ),
    class = "family"
  )
}

