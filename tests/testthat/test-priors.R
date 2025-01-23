
test_that("priors interface is working ", {
  #skip_on_cran()
  #skip_on_ci()

  # Requires loading RTMB
  library(RTMB)

  # Load data
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

  # Using fitRTMB
  log_prior = function(p){
    "c" <- ADoverload("c")
    "[<-" <- ADoverload("[<-")
    sum(dnorm( p$beta_z[9:16], mean=0, sd=0.25, log=TRUE))
  }
  fitRTMB = dsemRTMB( sem = sem,
              tsdata = Z,
              family = family,
              log_prior = log_prior,
              control = dsem_control( use_REML = FALSE) )

  # Run model using dsem ... not working in testthat mode
  neglog_prior = function(obj){
    "c" <- ADoverload("c")
    "[<-" <- ADoverload("[<-")
    -1 * sum(dnorm( obj$par[9:16], mean=0, sd=0.25, log=TRUE))
  }
  fit = dsem( sem=sem,
               tsdata=Z,
               family=family,
               prior_negloglike = neglog_prior,
               control = dsem_control(use_REML=FALSE) )

  # Check objective function
  expect_equal( as.numeric(fitRTMB$opt$obj), 198.1363, tolerance=1e-2 )
  expect_equal( as.numeric(fit$opt$obj), as.numeric(fitRTMB$opt$obj), tolerance=1e-2 )
})

