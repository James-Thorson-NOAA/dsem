

# R = (a S) / (1 + b S) * exp(e)
# 1/R = (1 + b S) / (a S) * exp(-e)
# S/R = (S + b S^2) / (a S) * exp(-e)
#     = ( (1/a) + (b/a) S ) * exp(-e)
# S/R ~ lognormal( meanlog, sdlog )
#   where meanlog = log( (1/a) + (b/a) S )
# i.e., `family = lognormal( link = "log" )`
#
# So
# slope: b / a = beta
# intercept:  1 / a = mu_1 + mu_2 * (b/a)

library(RTMB)

# Steepness parameterizations
R0 = 1
B0 = 1
h = 0.6
n = 100
sigmaR = 0.1

# Original parameters
a = 4*h*R0 / (1-h) / B0
b = (5*h - 1) / (1-h) / B0

# SR function
SR = function(B, pars, type = "steepness"){
  if(type=="steepness"){
    4 * pars$h * pars$R0 * B / (pars$B0 * (1-pars$h) + B * (5*pars$h-1))
  }else if(type=="a_b"){
    pars$a * B / ( 1 + pars$b * B )
  }
}

# Simulate values
B_t = runif( n, min = 0, max = 1.2 )
Rexp_t = SR( B_t, list(R0=R0, h=h, B0=B0) )
#Rexp_t = SR( B_t, list(a=a, b=b), type = "a_b" )
R_t = Rexp_t * exp( rnorm(n, sd = sigmaR) )

#####################
# Fit in RTMB
#####################

# Define nll
get_nll = function(p, data){
  getAll(p, data)
  Rhat = SR(B, pars = p, type="a_b")
  -sum(dnorm(log(R), log(Rhat), exp(logsigma),log=TRUE))
}

# Build obj
cmb <- function(f, d) function(p) f(p, d) ## Helper to make closure
obj = MakeADFun(
  func = cmb(
    get_nll,
    list( R = R_t, B = B_t )
  ),
  parameters = list(
    a = 1,
    b = 1,
    logsigma = 0
  )
)

# Fit
opt = nlminb( obj$par, obj$fn, obj$gr )

#####################
# Refit using dsem
#####################

#devtools::install_github("James-thorson-NOAA/dsem@dev_add-family")
library(dsem)
which_drop = sample(seq_len(n), size = 0.5*n, replace = FALSE)

# Define model
sem = "
  B -> BoverR, 0, BoverA, 1
  BoverR <-> BoverR, 0, NA, 0.001
  B <-> B, 0, sdB, 1
"

# Bundle data
df = data.frame(BoverR = B_t/R_t, B = B_t)
df$BoverR[which_drop] = NA

# Fit
fit = dsem(
  tsdata = ts(df),
  sem = sem,
  family = list( BoverR = lognormal("identity"), B = fixed() ),
  control = dsem_control(
    profile = NULL,
    newton_loops = 0,
    use_REML = FALSE,
    lower = c( 0.01, 0.01, log(0.001), 0.001, 0.001 ) # slope must be positive
  )
)

######################
# Compare parameter estimates
######################

parlist = as.list(fit$sdrep, what = "Estimate")
# Ratio
b / a
opt$par[2] / opt$par[1]
parlist$beta_z[1]
# Maximum
a
opt$par[1]
(ahat = 1 / (parlist$mu_j[1] - parlist$beta_z[1] * parlist$mu_j[2]))
# Density dependence
b
opt$par[2]
(bhat = parlist$beta_z[1] / (parlist$mu_j[1] - parlist$beta_z[1] * parlist$mu_j[2]))

# Plot
Rhat2 = SR( B, pars = list(a = ahat, b = bhat), type = "a_b")
matplot( x = B, y = cbind(Rhat, Rexp, Rhat2), type = "l", lty = "solid", lwd = 2, col = c("black","blue", "green") )
  points( x = B_t, y = R_t )


######################
# Or other dsem specification
######################

# Define model including intercept
sem = "
  B -> BoverR, 0, BoverA, 1
  ones -> BoverR, 0, b0
  BoverR <-> BoverR, 0, NA, 0
  B <-> B, 0, sdB, 1
  ones <-> ones, 0, NA, 1
"

# Bundle data
df = data.frame(BoverR = B_t/R_t, B = B_t, ones = 1)
df$BoverR[which_drop] = NA

# Fit
fit_b = dsem(
  tsdata = ts(df),
  sem = sem,
  family = list( BoverR = lognormal("identity"), B = fixed(), ones = fixed() ),
  estimate_mu = c("B"),
  control = dsem_control(
    profile = NULL,
    newton_loops = 0,
    #lower = c( 0.01, 0.01, log(0.001), 0.001, 0.001 ) # slope must be positive,
    use_REML = FALSE
  )
)

