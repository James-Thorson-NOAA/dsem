
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
  fit = dsem( 
    sem = sem,
    tsdata = tsdata,
    control = dsem_control(
      getJointPrecision = TRUE
    ) 
  )
  # Check objective function
  expect_equal( as.numeric(fit$opt$obj), 431.6122, tolerance=1e-2 )

  # Turning off RTMB for now
  if( FALSE ){
    fitRTMB = dsemRTMB( 
      sem = sem,
      tsdata = tsdata,
      control = dsem_control(
        getsd = FALSE
      ) 
    )
    # Check objective function
    expect_equal( as.numeric(fit$opt$obj), as.numeric(fitRTMB$opt$obj), tolerance=1e-2 )
  }
  
  # Convert and plot using phylopath
  as_fitted_DAG(fit)

  # Various other utilities
  plot(fit)
  vcov(fit, which="fixed")
  vcov(fit, which="random")
  vcov(fit, which="both")
  print(fit)
  logLik(fit)
  total_effect(fit)
  as_sem(fit)
  predict(fit, type="link")
  predict(fit, type="response")
  predict(fit, type="link", newdata=tsdata)
  simulate(fit, variance = "none")
  simulate(fit, variance = "random")
  simulate(fit, variance = "both")
  simulate(fit, resimulate_gmrf=TRUE)

  # Refit with measurement errors .. ignore quality of fit, just checking that it runs
  family = list(
    time = gaussian(),
    gnp = Gamma("log"),
    pwage = fixed(),
    cprofits = fixed(),
    consumption = fixed(),
    gwage = fixed(),
    invest = fixed(),
    capital = fixed()
  )
  fit1 = dsem(
    sem = sem,
    tsdata = tsdata,
    #family = c("normal","gamma",rep("fixed",ncol(tsdata)-2)),
    family = family,
    control = dsem_control(
      getsd=FALSE,
      extra_convergence_checks = FALSE,
      newton_loops=0
    )
  )
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
  fit_equations = dsem(
    sem = sem_equations,
    tsdata=tsdata,
    control = dsem_control(getsd=FALSE)
  )
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
  fit1 = dsem(
    sem = "",
    tsdata = data
  )
  # initial first without delta0 (to improve starting values)
  fit2 = dsem(
    sem = sem,
    tsdata = data
  )
  # Check objective function
  expect_equal( as.numeric(fit1$opt$obj), as.numeric(fit2$opt$obj), tolerance=1e-2 )
})

test_that("dsem works with fixed variances ", {
  data(isle_royale)
  data = ts( log(isle_royale[,2:3]), start=1959)

  # initial first without delta0 (to improve starting values)
  sem = "
    wolves <-> wolves, 0, sd1
    moose <-> moose, 0, sd2
    wolves -> wolves, 1, rho1
    moose -> moose, 1, rho2
  "
  fit1 = dsem(
    sem = sem,
    tsdata = data
  )

  # initial first without delta0 (to improve starting values)
  sem = "
    wolves <-> wolves, 0, NA, 0.3732812
    moose <-> moose, 0, NA, 0.1911209
    wolves -> wolves, 1, NA, 0.8558536
    moose -> moose, 1, NA, 0.9926087
  "
  fit2 = dsem(
    sem = sem,
    tsdata = data,
    family = list( wolves = gaussian(), moose = gaussian() )
  )
  expect_equal( as.numeric(fit1$opt$obj), as.numeric(fit2$opt$obj), tolerance=1e-2 )

  # initial first without delta0 (to improve starting values)
  # Turn off RTMB
  if( FALSE ){
    sem = "
      wolves <-> wolves, 0, NA, 0.3732812
      moose <-> moose, 0, NA, 0.1911209
      wolves -> wolves, 1, NA, 0.8558536
      moose -> moose, 1, NA, 0.9926087
    "
    fit3 = dsemRTMB( sem = sem,
                 tsdata = data,
                 family = c("normal", "normal") )
    expect_equal( as.numeric(fit1$opt$obj), as.numeric(fit3$opt$obj), tolerance=1e-2 )
  }
})

test_that("bering sea example is stable ", {
  data(bering_sea)

  #
  Z = ts( bering_sea )
  #family = rep( "fixed", ncol(bering_sea) )

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
  fit = dsem(
    sem = sem,
    tsdata = Z,
    #family = family,
    control = dsem_control(
      use_REML=FALSE
    )
  )

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
  fit = dsem(
    sem = sem,
    tsdata = tsdata,
    control = dsem_control(
      getsd=FALSE
    )
  )
  expect_equal( as.numeric(fit$opt$obj), 224.2993, tolerance=1e-2 )

  # Turning off RTMB
  if( FALSE ){
    fitRTMB = dsemRTMB( sem=sem,
                 tsdata=tsdata,
                 control = dsem_control(getsd=FALSE) )
    # Check objective function
    expect_equal( as.numeric(fit$opt$obj), as.numeric(fitRTMB$opt$obj), tolerance=1e-2 )
  }
})

test_that("dfa example is working ", {
  data( harborSealWA, package="MARSS")
  n_factors = 2

  # Add factors to data
  tsdata = harborSealWA[,c("SJI","EBays","SJF","PSnd","HC")]
  newcols = array( NA,
                   dim = c(nrow(tsdata),n_factors),
                   dimnames = list(NULL,paste0("F",seq_len(n_factors))) )
  tsdata = ts( cbind(tsdata, newcols), start=1978)
  
  # Scale and center (matches with MARSS does when fitting a DFA)
  tsdata = scale( tsdata, center=TRUE, scale=TRUE )
  
  # Automated version
  #sem = make_dfa( variables = c("SJI","EBays","SJF","PSnd","HC"),
  #                n_factors = n_factors )
  # Manual specification to show structure, using equations-and-lags interface
  equations = "
    # Loadings of variables onto factors
    SJI = L11(0.1) * F1
    EBays = L12(0.1) * F1 + L22(0.1) * F2
    SJF = L13(0.1) * F1 + L23(0.1) * F2
    PSnd = L14(0.1) * F1 + L24(0.1) * F2
    HC = L15(0.1) * F1 + L25(0.1) * F2
  
    # random walk for factors
    F1 = NA(1) * lag[F1,1]
    F2 = NA(1) * lag[F2,1]
  
    # Unit variance for factors
    V(F1) = NA(1)
    V(F2) = NA(1)
  
    # Zero residual variance for variables
    V(SJI) = NA(0)
    V(EBays) = NA(0)
    V(SJF) = NA(0)
    V(PSnd) = NA(0)
    V(HC) = NA(0)
  "
  sem = convert_equations(equations)
  
  # Initial fit
  family = Map(function(.) gaussian(), colnames(tsdata))
    family$F1 = fixed()
    family$F2 = fixed()
  mydsem0 = dsem(
    tsdata = tsdata,
    sem = sem,
    family = family,
    estimate_delta0 = TRUE,
    estimate_mu = vector(),
    control = dsem_control( 
      quiet = TRUE,
      gmrf_parameterization = "project",
      run_model = FALSE
    ) 
  )
  
  # fix all measurement errors at diagonal and equal
  map = mydsem0$tmb_inputs$map
  map$lnsigma_z = factor( rep(1,length(mydsem0$tmb_inputs$parameters$lnsigma_z)) )
  
  # Fix factors to have initial value, and variables to not
  map$delta0_j = factor( c(rep(NA,ncol(harborSealWA)-1), 1:n_factors) )
  
  # profile "delta0_j" to match MARSS (which treats initial condition as unpenalized random effect)
  mydfa = dsem( 
    tsdata = tsdata,
    sem = sem,
    family = family,
    estimate_delta0 = TRUE,
    estimate_mu = vector(),
    control = dsem_control( 
      quiet = TRUE,
      map = map,
      gmrf_parameterization = "project",
      use_REML = TRUE
    ) 
  )

  # Check objective function
  expect_equal( as.numeric(mydfa$opt$obj), 40.00618, tolerance=1e-2 )
})
