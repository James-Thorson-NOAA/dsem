

#devtools::install_github("James-thorson-NOAA/dsem@dev_add-family")
library(dsem)

n_t = 5
n_c = 3

sigF_t = seq( 0.1, 0.3, length = n_t )

eps_tc = matrix( rnorm(n_t*n_c), nrow = n_t, ncol = n_c )
eps_tc = sweep( eps_tc, MARGIN = 1, FUN = "*", STAT = sigF_t )

dat = data.frame(
  setNames( data.frame(eps_tc),letters[seq_len(n_c)]),
  "F" = NA,
  setNames( NA*data.frame(eps_tc),paste0("v_",letters[seq_len(n_c)]) )
)

sem = "
  a <-> a, 0, NA, 0
  b <-> b, 0, NA, 0
  c <-> c, 0, NA, 0
  #d <-> d, 0, NA, 0
  #e <-> e, 0, NA, 0
  v_a <-> v_a, 0, NA, 1
  v_b <-> v_b, 0, NA, 1
  v_c <-> v_c, 0, NA, 1
  #v_d <-> v_d, 0, NA, 1
  #v_e <-> v_e, 0, NA, 1
  F <-> F, 0, sdF, 0.1
  #F -> F, 1, NA, 1
  v_a -> a, 0, F
  v_b -> b, 0, F
  v_c -> c, 0, F
  #v_d -> d, 0, F
  #v_e -> e, 0, F
"

# Initial build
fit0 = dsem(
  tsdata = ts(dat),
  sem = sem,
  control = dsem_control(
    build_model = FALSE,
    use_REML = FALSE
  )
)

#
pars = fit0$parameters
pars$x_tj = ifelse( pars$x_tj==0, array(rnorm(prod(dim(pars$x_tj))),dim=dim(pars$x_tj)), pars$x_tj )

#
family = Map(function(.) fixed(), colnames(dat))
family$a = gaussian_fixed_sd( sd = 0.01 )
family$b = gaussian_fixed_sd( sd = 0.01 )
family$c = gaussian_fixed_sd( sd = 0.01 )
#family$d = gaussian_fixed_sd( sd = 0.01 )
#family$e = gaussian_fixed_sd( sd = 0.01 )

# Refit
fit = dsem(
  tsdata = ts(dat),
  sem = sem,
  family = family,
  estimate_mu = vector(),
  control = dsem_control(
    run_model = FALSE,
    use_REML = FALSE,
    parameters = pars,
    gmrf_parameterization = "project"
  )
)
fit$obj$fn( fit$obj$par )
fit$obj$gr( fit$obj$par )
rep = fit$obj$report()
