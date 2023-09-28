
#' Classify variables path
#'
#' \code{classify_variables} is copied from \code{sem:::classifyVariables}
#'
#' Copied with permission from John Fox under licence GPL (>= 2)
#'
#' @param model SEM model
#'
#' @return Tagged-list defining exogenous and endogenous variables
#' @export
classify_variables <-
function( model ){

    variables <- logical(0)
    for (paths in model[, 1]) {
        vars <- gsub(pattern=" ", replacement="", x=paths)
        vars <- sub("-*>", "->", sub("<-*", "<-", vars))
        if (grepl("<->", vars)) {
            vars <- strsplit(vars, "<->")[[1]]
            if (is.na(variables[vars[1]]))
                variables[vars[1]] <- FALSE
            if (is.na(variables[vars[2]]))
                variables[vars[2]] <- FALSE
        }
        else if (grepl("->", vars)) {
            vars <- strsplit(vars, "->")[[1]]
            if (is.na(variables[vars[1]]))
                variables[vars[1]] <- FALSE
            variables[vars[2]] <- TRUE
        }
        else if (grepl("<-", vars)) {
            vars <- strsplit(vars, "<-")[[1]]
            if (is.na(variables[vars[2]]))
                variables[vars[2]] <- FALSE
            variables[vars[1]] <- TRUE
        }
        else stop("incorrectly specified model")
    }
    list(endogenous = names(variables[variables]), exogenous = names(variables[!variables]))
}
