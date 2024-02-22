

test_that("dsem gmrf-parameterization options ", {
  data(isle_royale)
  data = ts( log(isle_royale[,2:3]), start=1959)
  
  sem = "
    # Link, lag, param_name
    wolves -> wolves, 1, arW
    moose -> wolves, 1, MtoW
    wolves -> moose, 1, WtoM
    moose -> moose, 1, arM
  "
  # initial build of object
  fit0 = dsem( sem = sem,
               tsdata = data,
               family = c("normal", "normal"),
               control = dsem_control(
                 nlminb_loops = 0,
                 newton_loops = 0,
                 getsd = FALSE) )
  #
  params = fit0$tmb_inputs$parameters
  params$lnsigma_j = log( c(0.1,0.1) )
  map = fit0$tmb_inputs$map
  map$lnsigma_j = factor( c(NA,NA) )
  
  # gmrf_parameterization = "separable"
  fit1 = dsem( sem = sem,
               tsdata = data,
               family = c("normal", "normal"),
               control = dsem_control(
                 nlminb_loops = 1,
                 newton_loops = 1,
                 getsd = TRUE,
                 map = map,
                 parameters = params,
                 gmrf_parameterization = "separable") )
  
  # gmrf_parameterization = "projection"
  fit2 = dsem( sem = sem,
               tsdata = data,
               family = c("normal", "normal"),
               control = dsem_control(
                 nlminb_loops = 1,
                 newton_loops = 1,
                 getsd = TRUE,
                 map = map,
                 parameters = fit1$obj$env$parList(),
                 gmrf_parameterization = "projection") )
  expect_equal( as.numeric(fit1$opt$obj), as.numeric(fit2$opt$obj), tolerance=1e-2 )
})

