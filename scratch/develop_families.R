
if( FALSE ){
  setwd(R'(C:\Users\James.Thorson\Desktop\Git\dsem\src\)')
  TMB::compile( "dsem.cpp" )
  setwd(R'(C:\Users\James.Thorson\Desktop\Git\dsem\)')
  devtools::document()
  devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\dsem)', force = TRUE )
}

set.seed(101)

# Simulate bovariate VAR
B = matrix( c(0.8, 0.1, -0.1, 0.8), 2, 2 )
x_ti = matrix( NA, nrow=100, ncol = 2)
x_ti[1,] = rnorm(n=2, mean = 0, sd = 1 )
for( t in 2:nrow(x_ti) ){
  x_ti[t,] = x_ti[t-1,] %*% B + rnorm(n=2, mean = 0, sd = 0.1)
}
x_ti = sweep( x_ti, MARGIN = 1, STATS = c(4,4), FUN = "+" )

# log-linked Poisson errors
y_ti = matrix( rpois( n = prod(dim(x_ti)),
               lambda = exp(x_ti) ), ncol = 2 )
# identity-linked normal errors
z_ti = matrix( rnorm( n = prod(dim(x_ti)),
               mean = x_ti, sd = 1 / exp(0.5*4) ), ncol = 2 )

#
sem = "
  x -> x, 1, b_xx
  x -> y, 1, b_xy
  y -> x, 1, b_yx
  y -> y, 1, b_yy
"
fit1 = dsem(
  tsdata = ts(data.frame(x=z_ti[,1], y=z_ti[,2])),
  sem = sem,
  family = list(x = fixed(), y = fixed() ),
  control = dsem_control(
    getsd = FALSE,
    newton_loops = 0
  )
)
