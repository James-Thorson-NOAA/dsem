% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{checkDepPackageVersion}
\alias{checkDepPackageVersion}
\title{Check for version mismatch in dependent binary packages}
\usage{
checkDepPackageVersion(
  dep_pkg = "TMB",
  this_pkg = "dsem",
  write_file = FALSE,
  warn = TRUE
)
}
\arguments{
\item{dep_pkg}{upstream package}

\item{this_pkg}{downstream package}

\item{write_file}{(logical) write version file and quit?}

\item{warn}{give warning?}
}
\value{
logical: TRUE if the binary versions match
}
\description{
Copied from glmmTMB with permission
}
