
test_that("dsem example is working ", {
  #skip_on_ci()
  data(KleinI, package="AER")
  TS = ts(data.frame(KleinI, "time"=time(KleinI) - 1931))

  # Specify by declaring each arrow and lag
  sem = "
    # Link, lag, param_name
    cprofits -> consumption, 0, a1
    cprofits -> consumption, 1, a2
    pwage -> consumption, 0, a3
    gwage -> consumption, 0, a3

    cprofits -> invest, 0, b1
    cprofits -> invest, 1, b2
    capital -> invest, 0, b3

    gnp -> pwage, 0, c2
    gnp -> pwage, 1, c3
    time -> pwage, 0, c1
  "
  tsdata = TS[,c("time","gnp","pwage","cprofits",'consumption',
                 "gwage","invest","capital")]

  # Fit model
  fit = dsem( sem=sem,
              tsdata=tsdata,
              control = dsem_control(getJointPrecision=TRUE) )
  # Check objective function
  expect_equal( as.numeric(fit$opt$obj), 431.6122, tolerance=1e-2 )

  #
  fitRTMB = dsemRTMB( sem=sem,
              tsdata=tsdata,
              control = dsem_control(getsd=FALSE) )
  # Check objective function
  expect_equal( as.numeric(fit$opt$obj), as.numeric(fitRTMB$opt$obj), tolerance=1e-2 )

  # Convert and plot using phylopath
  as_fitted_DAG(fit)

  # Various other utilities
  plot(fit)
  vcov(fit, which="fixed")
  vcov(fit, which="random")
  vcov(fit, which="both")
  print(fit)
  logLik(fit)
  as_sem(fit)
  predict(fit, type="link")
  predict(fit, type="response")
  predict(fit, type="link", newdata=tsdata)
  simulate(fit, variance = "none")
  simulate(fit, variance = "random")
  simulate(fit, variance = "both")
  simulate(fit, resimulate_gmrf=TRUE)

  # Refit with measurement errors .. ignore quality of fit, just checking that it runs
  fit1 = dsem( sem=sem,
               tsdata = tsdata,
               family = c("normal","gamma",rep("fixed",ncol(tsdata)-2)),
               control = dsem_control(getsd=FALSE,
                                      extra_convergence_checks = FALSE,
                                      newton_loops=0) )
  residuals(fit1, type="deviance")
  residuals(fit1, type="response")

  # Test equations including lags, starting values, and mirrored params
  equations = "
    consumption = a1(0.1)*cprofits + a2(0.1)*lag[cprofits,1]+ a3(0.1)*pwage + a3(0.2)*gwage
    invest = b1*cprofits + b2*lag[cprofits,1] + b3*capital
    pwage = c1*time + c2*gnp + c3*lag[gnp,1]
  "
  sem_equations = convert_equations(equations)
  # Fit model
  fit_equations = dsem( sem=sem_equations,
              tsdata=tsdata,
              control = dsem_control(getsd=FALSE) )
  # Check objective function
  expect_equal( as.numeric(fit$opt$obj), as.numeric(fit_equations$opt$obj), tolerance=1e-2 )
})

test_that("dsem adds variances ", {
  data(isle_royale)
  data = ts( log(isle_royale[,2:3]), start=1959)

  sem = "
    wolves <-> wolves, 0, sd1
    moose <-> moose, 0, sd2
  "
  # initial first without delta0 (to improve starting values)
  fit1 = dsem( sem = "",
               tsdata = data )
  # initial first without delta0 (to improve starting values)
  fit2 = dsem( sem = sem,
               tsdata = data )
  # Check objective function
  expect_equal( as.numeric(fit1$opt$obj), as.numeric(fit2$opt$obj), tolerance=1e-2 )
})

test_that("bering sea example is stable ", {
  data(bering_sea)

  #
  Z = ts( bering_sea )
  family = rep( "fixed", ncol(bering_sea) )

  # Specify model
  sem = "
    log_seaice -> log_CP, 0, seaice_to_CP
    log_CP -> log_Cfall, 0, CP_to_Cfall
    log_CP -> log_Esummer, 0, CP_to_E
    log_PercentEuph -> log_RperS, 0, Seuph_to_RperS
    log_PercentCop -> log_RperS, 0, Scop_to_RperS
    log_Esummer -> log_PercentEuph, 0, Esummer_to_Suph
    log_Cfall -> log_PercentCop, 0, Cfall_to_Scop
    SSB -> log_RperS, 0, SSB_to_RperS

    log_seaice -> log_seaice, 1, AR1, 0.001
    log_CP -> log_CP, 1,  AR2, 0.001
    log_Cfall -> log_Cfall, 1, AR4, 0.001
    log_Esummer -> log_Esummer, 1, AR5, 0.001
    SSB -> SSB, 1, AR6, 0.001
    log_RperS ->  log_RperS, 1, AR7, 0.001
    log_PercentEuph -> log_PercentEuph, 1, AR8, 0.001
    log_PercentCop -> log_PercentCop, 1, AR9, 0.001
  "

  # Run model
  fit = dsem( sem=sem,
               tsdata=Z,
               family=family,
               control = dsem_control(use_REML=FALSE) )

  # Check objective function
  expect_equal( as.numeric(fit$opt$obj), 189.3005, tolerance=1e-2 )
})

test_that("Fixing parameters works ", {
  sem = "
    A -> B, 0, NA, 1
    B -> C, 0, NA, 0.5
    # Extra links
    A -> D, 1, beta
    B -> D, 1, beta
  "
  set.seed(101)
  nobs = 40
  A = rnorm(nobs)
  B = A + rnorm(nobs)
  C = 0.5 * B + rnorm(nobs)
  D = rnorm(nobs)
  tsdata = ts(cbind(A=A, B=B, C=C, D=D))

  # Run models
  fit = dsem( sem=sem,
               tsdata=tsdata,
               control = dsem_control(getsd=FALSE) )
  fitRTMB = dsemRTMB( sem=sem,
               tsdata=tsdata,
               control = dsem_control(getsd=FALSE) )
  # Check objective function
  expect_equal( as.numeric(fit$opt$obj), 224.2993, tolerance=1e-2 )
  expect_equal( as.numeric(fit$opt$obj), as.numeric(fitRTMB$opt$obj), tolerance=1e-2 )

})