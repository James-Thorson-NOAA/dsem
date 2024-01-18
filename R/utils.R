#' @title Check for version mismatch in dependent binary packages
#' @description Copied from glmmTMB with permission
#' @param dep_pkg upstream package
#' @param this_pkg downstream package
#' @param write_file (logical) write version file and quit?
#' @param warn give warning?
#' @return logical: TRUE if the binary versions match
#' @importFrom utils packageVersion
#' @export
checkDepPackageVersion <- function(dep_pkg = "TMB",
                                   this_pkg = "dsem",
                                   write_file = FALSE,
                                   warn = TRUE) {
    cur_dep_version <- as.character(packageVersion(dep_pkg))
    fn <- sprintf("%s-version", dep_pkg)
    if (write_file) {
        cat(sprintf("current %s version=%s: writing file\n", dep_pkg, cur_dep_version))
        writeLines(cur_dep_version, con = fn)
        return(cur_dep_version)
    }
    fn <- system.file(fn, package=this_pkg)
    built_dep_version <- scan(file=fn, what=character(), quiet=TRUE)
    result_ok <- identical(built_dep_version, cur_dep_version)
    if(warn && !result_ok) {
        warning(
            "Package version inconsistency detected.\n",
            sprintf("%s was built with %s version %s",
                    this_pkg, dep_pkg, built_dep_version),
            "\n",
            sprintf("Current %s version is %s",
                    dep_pkg, cur_dep_version),
            "\n",
            sprintf("Please re-install %s from source ", this_pkg),
            "or restore original ",
            sQuote(dep_pkg), " package (see '?reinstalling' for more information)"
        )
    }
    return(result_ok)
}

#' @name reinstalling
#' @rdname reinstalling
#' @title Reinstalling binary dependencies
#'
#' @description The \code{dsem} package depends on several upstream packages, which it
#' uses in a way that depends heavily on their internal (binary) structure.
#' Sometimes, therefore, installing an update to one of these packages will
#' require that you re-install a \emph{binary-compatible} version of \code{dsem},
#' i.e. a version that has been compiled with the updated version of the upstream
#' package.
#' \itemize{
#' \item If you have development tools (compilers etc.) installed, you
#' should be able to re-install a binary-compatible version of the package by running
#' \code{install.packages("dsem", type="source")}. If you want to install
#' the development version of \code{dsem} instead, you can use
#' \code{remotes::install_github("James-Thorson-NOAA/dsem")}.
#' (On Windows, you can install development tools following the instructions at
#' \url{https://cran.r-project.org/bin/windows/Rtools/}; on MacOS, see
#' \url{https://mac.r-project.org/tools/}.)
#'
#' \item If you do \emph{not} have development tools and can't/don't want to
#' install them (and so can't install packages with compiled code from source),
#' you can revert the upstream package(s) to their previous binary version. For example, using the
#' \code{checkpoint} package:
#' \preformatted{
#' ## load (installing if necessary) the checkpoint package
#' while (!require("checkpoint")) install.packages("checkpoint")
#' ## retrieve build date of installed version of dsem
#' bd <- as.character(asDateBuilt(
#'       packageDescription("dsem",fields="Built")))
#' oldrepo <- getOption("repos")
#' use_mran_snapshot(bd) ## was setSnapshot() pre-checkpoint v1.0.0
#' install.packages("TMB")
#' options(repos=oldrepo) ## restore original repo
#' }
#' A similar recipe (substituting \code{Matrix} for \code{TMB} and \code{TMB} for \code{dsem})
#' can be used if you get warnings about an incompatibility between \code{TMB} and \code{Matrix}.
#' }
#' @details Copied from glmmTMB with permission
NULL

