
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
  fit = dsem( sem=sem, tsdata=Z )
  # Check objective function
  expect_equal( as.numeric(fit$opt$obj), 587.4755, tolerance=1e-2 )

  # Convert and plot using phylopath
  as_fitted_DAG(fit)

  # Various other utilities
  plot(fit)
  vcov(fit)
  print(fit)
  logLik(fit)
  as_sem(fit)
  predict(fit)
  predict(fit, newdata=Z)
})

