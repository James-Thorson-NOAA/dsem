###################
# Varying slopes
###################

n_t = 51
rho_x = 0.0
rho_y = 0
rho_b = 0.5
sd_x = 1
sd_y = 0.2
sd_b = 0.3
mu_b = 0.5

sim_ar1 = function(n, rho, sd){
  x = rnorm(n, sd = sd)
  for(t in 2:length(x)) x[t] = rho*x[t-1] + sqrt(1-rho^2)*x[t]
  return(x)
}

beta_t = mu_b + sim_ar1( n_t, rho_b, sd_b )
x_t = sim_ar1( n_t, rho_x, sd_x )
y_t = beta_t * x_t + sim_ar1( n_t, rho_y, sd_y )


if( FALSE ){
  library(KFAS)
  model <- SSModel(
    y ~ SSMtrend(degree = 1, Q = list(NA)) +
         SSMregression(~ x, Q = NA),
    data = data.frame(x = x_t, y = y_t )
  )
  fit <- fitSSM(model, inits = c(0.1, 0.1))
  coef(fit$model)
}

#######################
# AR1 intercept
#######################

# KFAS
trend <- function(sigma, n) {
  Z <- array(seq_len(n), c(1, 1, n))
  T <- R <- matrix(1, 1, 1)
  Q <- matrix(sigma^2, 1, 1)
  a1 <- 0
  P1 <- 10
  SSMcustom(Z, T, R, Q, a1, P1, n = n, state_names = "timevarying trend")
}

model <- SSModel(Nile ~ SSMbespoke(trend(NA, length(Nile))), H = NA)
updatefn <- function(pars, model){
  model$Q[1, 1, 1] <- exp(0.5 * pars[1])
  model$H[1, 1, 1] <- exp(0.5 * pars[2])
  model
}

fit <- fitSSM(model, c(1, 20), updatefn, method = "BFGS")
conf_intv <- predict(fit$model, interval = "confidence", level = 0.95)

ts.plot(
  cbind(Nile, conf_intv), 
  col = c(1, 2, 2, 2),
  ylab = "Predicted Annual flow", 
  main = "River Nile"
) 

# DSEM
library(dsem)
mydsem = dsem(
  tsdata = ts( data.frame(Nile = Nile), start = 1871),
  sem = "Nile -> Nile, 1, ar",
  family = "normal",
  control = dsem_control(
    use_REML = FALSE
  )
)

predict(mydsem)

matplot( 
  y = cbind(predict(mydsem), conf_intv[,1])
)

# MARSS
if( FALSE ){
  library(MARSS)
  nile = as.numeric(Nile)
  year = as.numeric(time(Nile))
  
  Z <- array(0, dim = c(1, 1, length(Nile)))
  Z[1, 1, ] <- year - mean(year)
  mod.list = list(Z = Z, A = matrix("a"), R = matrix("r"), B = matrix(1), 
      U = matrix(0), Q = matrix("q"), x0 = matrix("pi"))
  fit1 <- MARSS(nile, model = mod.list, silent = TRUE)
  fit2 <- MARSS(nile, model = mod.list, inits = fit1, method = "BFGS")
  ggplot2::autoplot(fit2, plot.type = "fitted.ytT")
}

x = x_t
y = y_t
  Z <- array(0, dim = c(1, 1, length(x)))
  Z[1, 1, ] <- x - mean(x)
  mod.list = list(Z = Z, A = matrix("a"), R = matrix("r"), B = matrix(1), 
      U = matrix(0), Q = matrix("q"), x0 = matrix("pi"))
  fit1 <- MARSS(y, model = mod.list, silent = TRUE)
  fit2 <- MARSS(y, model = mod.list, inits = fit1, method = "BFGS")
  ggplot2::autoplot(fit2, plot.type = "fitted.xtt1")
