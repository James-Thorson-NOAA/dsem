

#devtools::install_github("James-thorson-NOAA/dsem@dev_add-family")
#devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\dsem)' )
library(dsem)

n_t = 100
n_c = 3

sigF_t = seq( 0.1, 0.3, length = n_t )

eps_tc = matrix( rnorm(n_t*n_c), nrow = n_t, ncol = n_c )
eps_tc = sweep( eps_tc, MARGIN = 1, FUN = "*", STAT = sigF_t )

#########################
# Option-2
#########################


dat = data.frame(
  setNames( data.frame(eps_tc),letters[seq_len(n_c)]),
  F = NA,
  slope = scale( seq_len(n_t), center = TRUE, scale = TRUE )
)
dat$slope[ sample(seq_len(n_t),n_t/2) ] = NA

sem = "
  a <-> a, 0, F
  b <-> b, 0, F
  c <-> c, 0, F
  F <-> F, 0, sdF, 0.1
  F -> F, 1, NA, 1
  slope <-> slope, 0, sd_slope
  slope -> slope, 1, NA, 1
  slope -> F, 0, beta
"

# Initial build
fit = dsem(
  tsdata = ts(dat),
  sem = sem,
  control = dsem_control(
    build_model = TRUE,
    use_REML = FALSE,
    gmrf_parameterization = "full",
    logscale_moderating_variance = TRUE
  )
)

#
pars = fit0$parameters
pars$x_tj = ifelse( pars$x_tj==0, array(rnorm(prod(dim(pars$x_tj))),dim=dim(pars$x_tj)), pars$x_tj )
pars$x_tj[,(names(dat)=="F")] = abs(pars$x_tj[,(names(dat)=="F")])

# Using log-scale for moderating variances
#pars$mu_j[(names(dat)=="F")] = 1

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
  estimate_mu = c("F"),
  control = dsem_control(
    run_model = TRUE,
    use_REML = FALSE,
    parameters = pars,
    gmrf_parameterization = "full",
    logscale_moderating_variance = TRUE
  )
)
rep = fit$obj$report()


