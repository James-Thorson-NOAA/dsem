
test_that("dsem example is working ", {
  #skip_on_ci()
  sem = "
    Profits -> Consumption, 0, a2
    Profits -> Consumption, -1, a3
    Priv_wage -> Consumption, 0, a4
    Gov_wage -> Consumption, 0, a4
    Consumption <-> Consumption, 0, v1
    Consumption -> Consumption, -1, ar1
    Consumption -> Consumption, -2, ar2
    Profits -> Investment, 0, b2
    Profits -> Investment, -1, b3
    Capital_stock -> Investment, -1, b4
    Investment <-> Investment, 0, v2
    neg_Gov_wage <-> neg_Gov_wage, 0, v3
    GNP -> Priv_wage, 0, c2
    Taxes -> Priv_wage, 0, c2
    neg_Gov_wage -> Priv_wage, 0, c2
    GNP -> Priv_wage, -1, c3
    Taxes -> Priv_wage, -1, c3
    neg_Gov_wage -> Priv_wage, -1, c3
    Time -> Priv_wage, 0, c4
    Priv_wage <-> Priv_wage, 0, v4
    GNP <-> GNP, 0, v5
    Profits <-> Profits, 0, v6
    Capital_stock <-> Capital_stock, 0, v7
    Taxes <-> Taxes, 0, v8
    Time <-> Time, 0, v9
    Gov_wage <-> Gov_wage, 0, v10
    Gov_expense <-> Gov_expense, 0, v11
  "

  # Load data
  data(KleinI, package="AER")
  Data = as.data.frame(KleinI)
  Data = cbind( Data, "time" = seq(1,22)-11 )
  colnames(Data) = sapply( colnames(Data), FUN=switch,
             "consumption"="Consumption", "invest"="Investment",
             "cprofits"="Profits", "capital"="Capital_stock", "gwage"="Gov_wage",
             "pwage"="Priv_wage", "gexpenditure"="Gov_expense", "taxes"="Taxes",
             "time"="Time", "gnp"="GNP")
  Z = ts( cbind(Data, "neg_Gov_wage"=-1*Data[,'Gov_wage']) )

  # Fit model
  fit = dsem( sem=sem,
              tsdata=Z,
              control = dsem_control(getJointPrecision=TRUE) )
  # Check objective function
  expect_equal( as.numeric(fit$opt$obj), 587.4755, tolerance=1e-2 )

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
  predict(fit, type="link", newdata=Z)
  simulate(fit, variance = "none")
  simulate(fit, variance = "random")
  simulate(fit, variance = "both")
  simulate(fit, resimulate_gmrf=TRUE)

  # Refit with measurement errors
  fit1 = dsem( sem=sem,
               tsdata=Z,
               family = c("normal","gamma",rep("fixed",ncol(Z)-2)),
               control = dsem_control(getsd=FALSE, newton_loops=0) )
  residuals(fit1, type="deviance")
  residuals(fit1, type="response")
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

