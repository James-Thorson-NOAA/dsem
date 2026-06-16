

test_that("dsem gmrf-parameterization options ", {
  data(isle_royale)
  data = ts( log(isle_royale[,2:3]), start=1959)
  
  sem = "
    # Link, lag, param_name
    wolves -> wolves, 1, arW
    moose -> wolves, 1, MtoW
    wolves -> moose, 1, WtoM
    moose -> moose, 1, arM
    moose <-> wolves, 0, crosscor
  "
  # initial build of object
  fit0 = dsem( 
    sem = sem,
    tsdata = data,
    family = list( wolves = gaussian(), moose = gaussian() ),
    estimate_delta0 = TRUE,
    control = dsem_control(
      nlminb_loops = 0,
      newton_loops = 0,
      getsd = FALSE,
      extra = FALSE
    ) 
  )

  #
  params = fit0$tmb_inputs$parameters
  params$lnsigma_z = log( c(0.1,0.1) )
  map = fit0$tmb_inputs$map
  map$lnsigma_z = factor( c(NA,NA) )
  
  # gmrf_parameterization = "full"
  fit1 = dsem( 
    sem = sem,
    tsdata = data,
    family = list( wolves = gaussian(), moose = gaussian() ),
    estimate_delta0 = TRUE,
    control = dsem_control(
      map = map,
      parameters = params,
      nlminb_loops = 1,
      newton_loops = 0,
      getsd = TRUE,
      extra = FALSE,
      gmrf_parameterization = "full"
    ) 
  )
  
  # gmrf_parameterization = "projection"
  fit2 = dsem(
    sem = sem,
    tsdata = data,
    family = list( wolves = gaussian(), moose = gaussian() ),
    estimate_delta0 = TRUE,
    control = dsem_control(
      map = map,
      parameters = fit1$obj$env$parList(),
      gmrf_parameterization = "project"
    )
  )
  #expect_equal( as.numeric(fit1$opt$obj), as.numeric(fit2$opt$obj), tolerance=1e-2 )
  expect_equal( summary(fit1$sdrep), summary(fit1$sdrep), tolerance=1e-3 )
})

test_that("dsem works with high condition number ", {
  data(isle_royale)
  data = ts( log(isle_royale[,2:3]), start=1959)

  sem = "
    # Link, lag, param_name
    wolves -> wolves, 1, arW
    moose -> wolves, 1, MtoW
    wolves -> moose, 1, WtoM
    moose -> moose, 1, arM
    moose <-> wolves, 0, crosscor
  "
  # fit
  fit1 = dsem(
    sem = sem,
    tsdata = data,
    #family = c("normal", "normal"),
    estimate_delta0 = TRUE,
    control = dsem_control(
      gmrf_parameterization = "gmrf_project",
      newton_loops = 0,
      getsd = FALSE,
      extra = FALSE
    )
  )
  # fit
  fit2 = dsem(
    sem = sem,
    tsdata = data,
    #family = c("normal", "normal"),
    estimate_delta0 = TRUE,
    control = dsem_control(
      gmrf_parameterization = "full",
      newton_loops = 0,
      getsd = FALSE,
      extra = FALSE
    )
  )
  expect_equal( as.numeric(fit1$opt$obj), 5.020765, tolerance=1e-3 )
  expect_equal( as.numeric(fit2$opt$obj), 5.020765, tolerance=1e-3 )
})

test_that("dsem constant-variance options ", {
  data(isle_royale)
  data = ts( log(isle_royale[,2:3]), start=1959)
  
  # Show that constant_variance = "marginal" has constant marginal variance *with* crosscorrelation
  sem = "
    # Link, lag, param_name
    wolves -> wolves, 1, arW
    moose -> wolves, 1, MtoW
    wolves -> moose, 1, WtoM
    moose -> moose, 1, arM
    moose <-> wolves, 0, NA, 0.2
  "
  # initial build of object
  fit = dsem( 
    sem = sem,
    tsdata = data,
    estimate_delta0 = TRUE,
    control = dsem_control(
      nlminb_loops = 1,
      newton_loops = 0,
      getsd = FALSE,
      constant_variance = "marginal"
    ) 
  )
  margvar = array( diag(as.matrix(solve(fit$obj$report()$Q_oo))), dim=dim(data))
  expect_equal( apply(margvar,MARGIN=2,FUN=sd), c(0,0), tolerance=0.05 )

  # Show that constant_variance = "diagonal" has constant marginal variance *without* crosscorrelation 
  sem = "
    # Link, lag, param_name
    wolves -> wolves, 1, arW
    moose -> wolves, 1, MtoW
    wolves -> moose, 1, WtoM
    moose -> moose, 1, arM
  "
  fit = dsem( 
    sem = sem,
    tsdata = data,
    estimate_delta0 = TRUE,
    control = dsem_control(
      nlminb_loops = 1,
      newton_loops = 0,
      getsd = FALSE,
      constant_variance = "diagonal"
    ) 
  )
  margvar = array( diag(as.matrix(solve(fit$obj$report()$Q_oo))), dim=dim(data))
  expect_equal( apply(margvar,MARGIN=2,FUN=sd), c(0,0), tolerance=0.01 )

  # Show that marginal variance increases
  sem = "
    # Link, lag, param_name
    wolves -> wolves, 1, arW
    moose -> wolves, 1, MtoW
    wolves -> moose, 1, WtoM
    moose -> moose, 1, arM
  "
  fit0 = dsem( 
    sem = sem,
    tsdata = data,
    estimate_delta0 = FALSE,
    control = dsem_control(
      nlminb_loops = 1,
      newton_loops = 0,
      getsd = FALSE,
      constant_variance = "conditional"
    ) 
  )
  parameters = fit0$obj$env$parList()
    parameters$delta0_j = rep( 0, ncol(data) )
  fit = dsem( 
    sem = sem,
    tsdata = data,
    estimate_delta0 = TRUE,
    control = dsem_control(
      nlminb_loops = 1,
      newton_loops = 0,
      getsd = FALSE,
      constant_variance = "conditional",
      parameters = parameters
    ) 
  )
  margvar = array( diag(as.matrix(solve(fit$obj$report()$Q_oo))), dim=dim(data))
})


#test_that("dfa using dsem is working ", {
#  data(isle_royale)
#  data = ts( cbind(log(isle_royale[,2:3]), "F"=NA), start=1959)
#  
#  sem = "
#    F -> wolves, 0, l1
#    F -> moose, 0, l2
#    F -> F, 1, NA, 1
#    F <-> F, 0, NA, 1
#    wolves <-> wolves, 0, NA, 0.01
#    moose <-> moose, 0, NA, 0.01
#  "
#  # initial build of object
#  fit = dsem( sem = sem,
#               tsdata = data,
#               family = c("normal", "normal", "normal"),
#               estimate_delta0 = FALSE,
#               control = dsem_control(
#                 gmrf_parameterization = "full",
#                 run_model = FALSE) )
#  Report = fit$obj$report()
#  
#  library(Matrix)
#  image(Matrix(solve(Report$IminusRho_kk)))
#})

test_that("dsem `gmrf_project` and `mvn_project` are working ", {

  make_ar = function(rho, X){
    for(t in 2:length(X)) X[t] = rho * X[t-1] + sqrt(1-rho) * X[t]
    return(X)
  }
  
  set.seed(123)
  X = rnorm(100)
  X = make_ar( rho = 0.8, X = X )
  p = plogis(X)
  Y = rbinom( n = length(p), size = 1, prob = p )
  
  # Bundle
  dat = data.frame( 
    X = X, 
    Y = Y 
  )
  
  sem = "
    X -> X, 1, rho
    X -> Y, 0, b_XY
    Y <-> Y, 0, NA, 0
    Y -> X, 2, gamma     # Not included but forces loops for testing
  "
  # New option
  control = dsem_control(
    gmrf_parameterization = "gmrf_project",
    use_REML = FALSE,
    newton_loops = 0
  )
  fit1 = dsem(
    tsdata = ts(dat),
    sem = sem,
    control = control,
    family = list( X = fixed(), Y = binomial("logit") )
  )
  
  # New option
  control = dsem_control(
    gmrf_parameterization = "mvn_project",
    use_REML = FALSE
  )
  fit2 = dsem(
    tsdata = ts(dat),
    sem = sem,
    control = control,
    family = list( X = fixed(), Y = binomial("logit") )
  )

  # Old option
  sem = "
    X -> X, 1, rho
    X -> Y, 0, b_XY
    Y <-> Y, 0, NA, 0.0001
    Y -> X, 2, gamma
  "
  control = dsem_control(
    use_REML = FALSE,
    #newton_loops = 0,
    #nlminb_loops = 1,
    #getsd = FALSE,
    #extra = FALSE,
    gmrf_parameterization = "full"
  )
  fit0 = dsem(
    tsdata = ts(dat),
    sem = sem,
    control = control,
    family = list( X = fixed(), Y = binomial("logit") )
  )

  # `gmrf_project` without any projection
  control = dsem_control(
    gmrf_parameterization = "gmrf_project",
    use_REML = FALSE,
    newton_loops = 0
  )
  fit3 = dsem(
    tsdata = ts(dat),
    sem = sem,
    control = control,
    family = list( X = fixed(), Y = binomial("logit") )
  )

  expect_equal( summary(fit1), summary(fit2), tolerance=0.001 )
  expect_equal( summary(fit1), summary(fit0), tolerance=0.001 )
  expect_equal( summary(fit1), summary(fit3), tolerance=0.001 )
})
