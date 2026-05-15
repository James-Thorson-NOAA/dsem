
test_that("random slopes example is working ", {
  # Load data
  data(pdo_departure_bay)
  
  # Format
  tsdata = ts(data.frame(
    Temp = pdo_departure_bay[,2],
    PDO = pdo_departure_bay[,3],
    slope = NA
  ), start = 1914 )
  
  # Model
  sem = "
    PDO -> Temp, 0, slope
    slope -> slope, 1, ar_slope
    PDO -> PDO, 1, ar_PDO
    Temp -> Temp, 1, ar_Temp
  "
  
  # Fit
  fit = dsem(
    tsdata = tsdata,
    sem = sem,
    estimate_mu = colnames(tsdata),
    estimate_delta0 = FALSE,
    control = dsem_control(
      quiet = TRUE,
      use_REML = FALSE,
      gmrf_parameterization = "full"
    )
  )

  # Check objective function
  expect_equal( as.numeric(fit$opt$obj), 245.6308, tolerance=1e-2 )
})

test_that("Lotka-Volterra example is working ", {
  data(paramesium_didinium)
  orig_dat = data.frame( 
    X = paramesium_didinium[,'paramecium'] / 100,
    Y = paramesium_didinium[,'didinium'] / 100 
  )

  # Format
  dat = full_dat = cbind(
    logX = log(orig_dat$X), logY = log(orig_dat$Y),
    X = NA, Y = NA,
    logX1 = NA, logY1 = NA,
    logX2 = NA, logY2 = NA,
    logX3 = NA, logY3 = NA,
    ones = 1
  )
  
  # Center variables for numerical stability
  mean_j = colMeans( dat[,1:2], na.rm = TRUE )
  dat[,1:2] = sweep( dat[,1:2], FUN = "-", MARGIN = 2, STATS = mean_j )
  
  sem = "
    # Main interactions
    logX -> logX, 1, NA, 1
    ones -> logX, 0, alpha
    Y -> logX, 1, beta, -0.1
  
    # Form X \approx exp(logX)
    ones -> X, 0, NA, 1
    logX -> logX1, 0, NA, 1
    logX1 -> X, 0, NA, 1
    logX1 -> logX2, 0, logX
    logX2 -> X, 0, NA, 0.5
    logX2 -> logX3, 0, logX
    logX3 -> X, 0, NA, 0.166
  
    # Variances
    X <-> X, 0, NA, 0
    logX <-> logX, 0, sd_logX
    logX1 <-> logX1, 0, NA, 0
    logX2 <-> logX2, 0, NA, 0
    logX3 <-> logX3, 0, NA, 0
  
    # Main interactions
    logY -> logY, 1, NA, 1
    X -> logY, 1, gamma
    ones -> logY, 0, delta, -0.1
  
    # Form Y \approx exp(logY)
    ones -> Y, 0, NA, 1
    logY -> logY1, 0, NA, 1
    logY1 -> Y, 0, NA, 1
    logY1 -> logY2, 0, logY
    logY2 -> Y, 0, NA, 0.5
    logY2 -> logY3, 0, logY
    logY3 -> Y, 0, NA, 0.166
  
    # Variances
    Y <-> Y, 0, NA, 0
    logY <-> logY, 0, sd_logY
    logY1 <-> logY1, 0, NA, 0
    logY2 <-> logY2, 0, NA, 0
    logY3 <-> logY3, 0, NA, 0
  
    # Dummy constant
    ones <-> ones, 0, NA, 0.001
    ones -> ones, 1, NA, 1
  "
  
  fit = dsem(
    tsdata = ts(dat),
    sem = sem,
    estimate_mu = vector(), 
    estimate_delta0 = FALSE,
    control = dsem_control(
      quiet = TRUE
    )
  )

  # Check objective function
  expect_equal( as.numeric(fit$opt$obj), 499647.8008, tolerance=1e-2 )
})
