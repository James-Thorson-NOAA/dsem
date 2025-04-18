% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_equations.R
\name{convert_equations}
\alias{convert_equations}
\title{Convert equations notation}
\usage{
convert_equations(equations)
}
\arguments{
\item{equations}{Specification for time-series structural equation model structure
including lagged or simultaneous effects.  See Details section in
\code{\link[dsem]{convert_equations}} for more description}
}
\description{
Converts equations to arrow-and-lag notation expected by dsem
}
\details{
The function modifies code copied from package
`sem` under licence GPL (>= 2) with permission from John Fox.

For specifyEquations, each input line is either a regression equation or the
specification of a variance or covariance. Regression equations are of the form
y = par1*x1 + par2*x2 + ... + park*xk
where y and the xs are variables in the model (either observed or latent),
and the pars are parameters. If a parameter is given as a numeric value
(e.g., 1) then it is treated as fixed. Note that no error variable is
included in the equation; error variances are specified via either
the covs argument, via V(y) = par (see immediately below), or are
added automatically to the model when, as by default, endog.variances=TRUE.
A regression equation may be split over more than one input by breaking at a +,
so that + is either the last non-blank character on a line or the
first non-blank character on the subsequent line.

Variances are specified in the form V(var) = par and
covariances in the form C(var1, var2) = par, where the vars are
variables (observed or unobserved) in the model. The symbols V and C
may be in either lower- or upper-case. If par is a numeric value (e.g., 1)
then it is treated as fixed. In conformity with the RAM model,
a variance or covariance for an endogenous variable in the
model is an error variance or covariance.

To set a start value for a free parameter, enclose the numeric
start value in parentheses after the parameter name, as parameter(value).
}
