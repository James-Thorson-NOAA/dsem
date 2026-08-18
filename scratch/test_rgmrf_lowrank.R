


# Simulate a log-linked Poisson GLM
x = rnorm(100)
logodds = 1 + 0.8 * x
p = plogis(logodds)
y = rbinom( length(x), prob = p, size = 1 )
data = data.frame(x=x, y=y)

# Define zero process errors (only process errors) in response y
sem = "
  x <-> y, 0, beta
  y <-> y, 0, NA, 0
"

# Fit as DSEM
fit = dsem(
  sem = sem,
  tsdata = ts(data),
  family = list( x = fixed(), y = binomial("logit") ),
  estimate_mu = c("y"),
  control = dsem_control(quiet=TRUE, use_REML = FALSE)
)


simulate(fit, resimulate_gmrf = TRUE)

# Extract covariance
M = fit$obj$report()$IminusRho_kk
G = fit$obj$report()$Gamma
Q = t(M) %*% solve(t(G) %*% G) %*% M
Sigma1 = solve(as.matrix(Q))

##
#x_kz = rgmrf_lowrank(
#  G = G,
#  M = M,
#  mu = rep(0,nrow(M))
#)
#cov(t(as.matrix(x_kz)))

