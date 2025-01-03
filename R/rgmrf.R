rgmrf <-
function( mu, # estimated fixed and random effects
          Q # estimated joint precision
          ) {

  # Simulate values
  if(missing(mu)) mu = rep(0,nrow(Q))
  z0 = rnorm(length(mu))
  # Q = t(P) * L * t(L) * P
  L = Matrix::Cholesky(Q, super=TRUE)
  # Calcualte t(P) * solve(t(L)) * z0 in two steps
  z = Matrix::solve(L, z0, system = "Lt") # z = Lt^-1 * z
  z = Matrix::solve(L, z, system = "Pt") # z = Pt    * z
  return(mu + z)
}
