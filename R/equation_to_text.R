#' @title Convert equations notation
#'
#' @description Converts equations to arrow-and-lag notation expected by dsem
#'
#' @param equations Specification for time-series structural equation model structure
#'        including lagged or simultaneous effects.  See Details section in
#'        \code{\link[dsem]{equations_to_text}} for more description
#'
#' @details
#' The function modifies code copied from package
#' `sem` under licence GPL (>= 2) with permission from John Fox.
#'
#' For specifyEquations, each input line is either a regression equation or the
#' specification of a variance or covariance. Regression equations are of the form
#' y = par1*x1 + par2*x2 + ... + park*xk
#' where y and the xs are variables in the model (either observed or latent),
#' and the pars are parameters. If a parameter is given as a numeric value
#' (e.g., 1) then it is treated as fixed. Note that no error variable is
#' included in the equation; error variances are specified via either
#' the covs argument, via V(y) = par (see immediately below), or are
#' added automatically to the model when, as by default, endog.variances=TRUE.
#' A regression equation may be split over more than one input by breaking at a +,
#' so that + is either the last non-blank character on a line or the
#' first non-blank character on the subsequent line.
#'
#' Variances are specified in the form V(var) = par and
#' covariances in the form C(var1, var2) = par, where the vars are
#' variables (observed or unobserved) in the model. The symbols V and C
#' may be in either lower- or upper-case. If par is a numeric value (e.g., 1)
#' then it is treated as fixed. In conformity with the RAM model,
#' a variance or covariance for an endogenous variable in the
#' model is an error variance or covariance.
#'
#' To set a start value for a free parameter, enclose the numeric
#' start value in parentheses after the parameter name, as parameter(value).
#'
#'
#' @export
equation_to_text <-
function(equations){

  # Local functions
  not.number <- function (constant){
    save <- options(warn = -1)
    on.exit(save)
    is.na(as.numeric(constant))
  }
  par.start <- function(coef, eq) {
    if (length(grep("\\(", coef)) == 0) {
        return(c(coef, "NA"))
    }
    par.start <- strsplit(coef, "\\(")[[1]]
    if (length(par.start) != 2)
        stop("Parse error in equation: ", eq, "\n  Start values must be given in the form \"parameter(value)\".")
    par <- par.start[[1]]
    start <- par.start[[2]]
    if (length(grep("\\)$", start)) == 0)
        stop("Parse error in equation: ", eq, "\n  Unbalanced parentheses.")
    start <- sub("\\)", "", start)
    return(c(par, start))
  }
  parseEquation <- function(eqn) {
    eq <- eqn
    eqn <- gsub("\\s*", "", eqn)
    eqn <- strsplit(eqn, "=")[[1]]
    if (length(eqn) != 2)
        stop("Parse error in equation: ", eq, "\n  An equation must have a left- and right-hand side separated by =.")
    lhs <- eqn[1]
    rhs <- eqn[2]
    if (length(grep("^[cC]\\(", lhs)) > 0) {
        if (length(grep("\\)$", lhs)) == 0)
            stop("Parse error in equation: ", eq, "\n  Unbalanced parentheses.")
        lhs <- sub("[cC]\\(", "", lhs)
        lhs <- sub("\\)", "", lhs)
        variables <- strsplit(lhs, ",")[[1]]
        if (length(variables) != 2)
            stop("Parse error in equation: ", eq, "\n  A covariance must be in the form C(var1, var2) = cov12")
        if (not.number(rhs)) {
            par.start <- par.start(rhs, eq)
            if (not.number(par.start[2]) && (par.start[2] !=
              "NA"))
              stop("Parse error in equation: ", eq, "\n  Start values must be numeric constants.")
            ram <- paste(variables[1], " <-> ", variables[2],
              ", ", par.start[1], ", ", par.start[2], sep = "")
        }
        else {
            ram <- paste(variables[1], " <-> ", variables[2],
              ", NA, ", rhs, sep = "")
        }
    }
    else if (length(grep("^[vV]\\(", lhs)) > 0) {
        lhs <- sub("[vV]\\(", "", lhs)
        if (length(grep("\\)$", lhs)) == 0)
            stop("Parse error in equation: ", eq, "\n  Unbalanced parentheses.")
        lhs <- sub("\\)", "", lhs)
        if (not.number(rhs)) {
            par.start <- par.start(rhs, eq)
            if (not.number(par.start[2]) && (par.start[2] !=
              "NA"))
              stop("Parse error in equation: ", eq, "\n  Start values must be numeric constants.")
            ram <- paste(lhs, " <-> ", lhs, ", ", par.start[1],
              ", ", par.start[2], sep = "")
        }
        else {
            ram <- paste(lhs, " <-> ", lhs, ", NA, ", rhs,
              sep = "")
        }
    }
    else {
        terms <- strsplit(rhs, "\\+")[[1]]
        terms <- strsplit(terms, "\\*")
        ram <- character(length(terms))
        for (term in 1:length(terms)) {
            trm <- terms[[term]]
            if (length(trm) != 2)
              stop("Parse error in equation: ", eq, "\n  The term  \"",
                trm, "\" is malformed.", "\n  Each term on the right-hand side of a structural equation must be of the form \"parameter*variable\".")
            coef <- trm[1]
            if (not.number(coef)) {
              par.start <- par.start(coef, eq)
              if (not.number(par.start[2]) && (par.start[2] !=
                "NA"))
                stop("Parse error in equation: ", eq, "\n  Start values must be numeric constants.")
              ram[term] <- paste(trm[2], " -> ", lhs, ", ",
                par.start[1], ", ", par.start[2], sep = "")
            }
            else {
              ram[term] <- paste(trm[2], " -> ", lhs, ", NA, ",
                coef, sep = "")
            }
        }
    }
    ram
  }
  parse_lag <- function(term){
    term_split = strsplit( term, " -> ", fixed=TRUE )[[1]]
    var = term_split[1]
    if (length(grep("\\[", var)) == 0) {
      return(c(term, "0"))
    }
    var_split <- strsplit(var, "\\[")[[1]]
    if (length(var_split) != 2){
      stop("Parse error in equation: ", term, "\n  Lags must be given in the form \"lag[var,integer]\".")
    }
    par_split = strsplit(var_split[[2]], ",", fixed=TRUE)[[1]]
    par = par_split[1]
    lag <- par_split[[2]]
    if (length(grep("\\]$", lag)) == 0){
      stop("Parse error in equation: ", term, "\n  Unbalanced parentheses.")
    }
    lag <- sub("\\]", "", lag)
    term_out = paste0( par, " -> ", term_split[2] )
    return(c(term_out, lag))
  }
  add_lags <- function(term){
    term_split = strsplit(term, ", ", fixed=TRUE)[[1]]
    lag = parse_lag(term_split[1])
    term_ram = paste0( lag[1], ", ", lag[2], ", ", term_split[2], ", ", term_split[3] )
    return(term_ram)
  }

  # Read text
  equations <- scan(text = equations, what = "", sep = ";", strip.white = TRUE, comment.char = "#")
  equations2 <- character(0)
  eqn <- 0
  skip <- FALSE
  for (equation in equations) {
      eqn <- eqn + 1
      if (skip) {
          skip <- FALSE
          next
      }
      if (substring(equation, 1, 1) == "+") {
          equations2[length(equations2)] <- paste(equations2[length(equations2)],
              equation)
      }
      else if (substring(equation, nchar(equation)) == "+") {
          equations2 <- c(equations2, paste(equation, equations[eqn +
              1]))
          skip <- TRUE
      }
      else equations2 <- c(equations2, equation)
  }
  ram <- unlist(lapply(equations2, parseEquation))

  #
  ram <- unlist(lapply(ram, add_lags))
  return(ram)
}
