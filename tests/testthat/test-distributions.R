
test_that("dsem example is working ", {
  #devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\dsem)', force=TRUE, dep=FALSE )  
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
  fit1 = dsem( tsdata = ts(data.frame(x=z_ti[,1], y=z_ti[,2])),
        sem = sem,
        family = c("fixed", "fixed") )
  fit2 = dsem( tsdata = ts(data.frame(x=y_ti[,1], y=y_ti[,2])),
        sem = sem,
        family = c("poisson", "poisson") )
  fit3 = dsemRTMB( tsdata = ts(data.frame(x=y_ti[,1], y=y_ti[,2])),
        sem = sem,
        family = c("poisson", "poisson") ) 
  expect_equal( as.numeric(fit2$opt$obj), as.numeric(fit2$opt$obj), tolerance=1e-2 )
  expect_equal( resid(fit2), resid(fit3), tolerance=1e-2 )
})

