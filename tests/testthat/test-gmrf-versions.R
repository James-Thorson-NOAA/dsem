

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
  fit0 = dsem( sem = sem,
               tsdata = data,
               family = c("normal", "normal"),
               estimate_delta0 = TRUE,
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
               estimate_delta0 = TRUE,
               control = dsem_control(
                 map = map,
                 parameters = params,
                 gmrf_parameterization = "separable") )
  
  # gmrf_parameterization = "projection"
  fit2 = dsem( sem = sem,
               tsdata = data,
               family = c("normal", "normal"),
               estimate_delta0 = TRUE,
               control = dsem_control(
                 map = map,
                 parameters = fit1$obj$env$parList(),
                 gmrf_parameterization = "projection") )
  expect_equal( as.numeric(fit1$opt$obj), as.numeric(fit2$opt$obj), tolerance=1e-2 )
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
#                 gmrf_parameterization = "separable",
#                 run_model = FALSE) )
#  Report = fit$obj$report()
#  
#  library(Matrix)
#  image(Matrix(solve(Report$IminusRho_kk)))
#})
