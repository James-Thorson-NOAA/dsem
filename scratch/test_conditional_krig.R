
# devtools::install_local( R'(C:\Users\james\OneDrive\Desktop\Git\dsem)', force = TRUE, dep = FALSE )

X = cumsum( rnorm(100) )
Y = 2 * X
p_missing = 0.2

Xobs = X
Xobs[sample(seq_along(X),size=p_missing*length(X),replace=FALSE)] = NA

dat = data.frame( X = Xobs, Y = NA )

library(dsem)

sem = "
  X -> X, 1, rho
  X -> Y, 0, NA, 2
  Y <-> Y, 0, NA, 0
"

control = dsem_control(
  run_model = TRUE,
  use_REML = FALSE
)
control$gmrf_parameterization = "mvn_project"
fit = dsem(
  tsdata = ts(dat),
  sem = sem,
  control = control
)
rep = fit$obj$report()

obj = fit$obj
