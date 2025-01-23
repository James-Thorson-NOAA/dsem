
library(RTMB)

# Works
func <- function(p) {
  u <- p$u
  ans <- -dnorm(u[1], log=TRUE) ## u[1] ~ N(0,1)
  ans <- ans - sum(dnorm(diff(u), log=TRUE)) ## u[i]-u[i-1] ~ N(0,1)
  x = outer(u, rep(1,10))
  REPORT( x )
  return(ans)
}
obj <- MakeADFun(func, list(u=numeric(20)), random="u")
obj$simulate()

# Works
func <- function(p) {
  u <- p$u
  x = u
  REPORT( x )
  ans <- -dnorm(u[1], log=TRUE) ## u[1] ~ N(0,1)
  ans <- ans - sum(dnorm(diff(u), log=TRUE)) ## u[i]-u[i-1] ~ N(0,1)
  return(ans)
}
obj <- MakeADFun(func, list(u=numeric(20)), random="u")
obj$simulate()

# Fails
func <- function(p) {
  u <- p$u
  x = outer(u, rep(1,2))
  REPORT( x )
  ans <- -dnorm(u[1], log=TRUE) ## u[1] ~ N(0,1)
  ans <- ans - sum(dnorm(diff(u), log=TRUE)) ## u[i]-u[i-1] ~ N(0,1)
  return(ans)
}
obj <- MakeADFun(func, list(u=numeric(20)), random="u")
obj$simulate()

# Fails
func <- function(p) {
  u <- p$u
  x = matrix( c(u,u), ncol=2 )
  REPORT( x )
  ans <- -dnorm(u[1], log=TRUE) ## u[1] ~ N(0,1)
  ans <- ans - sum(dnorm(diff(u), log=TRUE)) ## u[i]-u[i-1] ~ N(0,1)
  return(ans)
}
obj <- MakeADFun(func, list(u=numeric(20)), random="u")
obj$simulate()



