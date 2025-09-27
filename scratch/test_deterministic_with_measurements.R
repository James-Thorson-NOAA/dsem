
# devtools::install_local( R'(C:\Users\james\OneDrive\Desktop\Git\dsem)', force = TRUE, dep = FALSE )

library(dsem)

make_ar = function(rho, X){
  for(t in 2:length(X)) X[t] = rho * X[t-1] + sqrt(1-rho) * X[t]
  return(X)
}

################
# 
# USE-CASE
#  Determinstic and family = "fixed" with measurements
#
################


p_missing = 0.4

X = cumsum( rnorm(100) )
Y = 2 * X
Z = 0.5*Y + rnorm(100)
X_missing = sample(seq_along(X), size=p_missing*length(X), replace=FALSE)
Y_missing = sample(seq_along(X), size=p_missing*length(X), replace=FALSE)
Y_missing = union( setdiff( seq_along(X), X_missing ), Y_missing )

# Missing-ness
Xobs = X
Xobs[X_missing] = NA
Yobs = Y
Yobs[Y_missing] = NA

# Bundle
dat = data.frame( 
  X = 4 + Xobs, 
  Y = 1 + Yobs, 
  Z = 2 + Z 
)

library(dsem)

sem = "
  X -> X, 1, rho
  X -> Y, 0, b_XY
  Y -> Z, 0, b_YZ
  Y <-> Y, 0, NA, 0
"

control = dsem_control(
  gmrf_parameterization = "mvn_project",
  project_k = is.na(dat)
)
fit = dsem(
  tsdata = ts(dat),
  sem = sem,
  control = control
)

################
# 
# USE-CASE
#  logistic regression
#
################

X = rnorm(100)
X = make_ar( rho = 0.8, X = X )
p = plogis(X)
Y = rbinom( n = length(p), size = 1, prob = p )

# Bundle
dat = data.frame( 
  X = X, 
  Y = Y 
)


sem = "
  X -> X, 1, rho
  X -> Y, 0, b_XY
  Y <-> Y, 0, NA, 0
"

# New option
control = dsem_control(
  gmrf_parameterization = "mvn_project",
  #build_model = FALSE,
  use_REML = FALSE
)
fit = dsem(
  tsdata = ts(dat),
  sem = sem,
  control = control,
  family = c("fixed", "bernoulli")
)

# Old option
sem = "
  X -> X, 1, rho
  X -> Y, 0, b_XY
  Y <-> Y, 0, NA, 0.0001
"
control = dsem_control(
  use_REML = FALSE
)
fit0 = dsem(
  tsdata = ts(dat),
  sem = sem,
  control = control,
  family = c("fixed", "bernoulli")
)

# Other new option
sem = "
  X -> X, 1, rho
  X -> Y, 0, b_XY
  Y <-> Y, 0, NA, 0
"
control = dsem_control(
  use_REML = FALSE,
  stabilize_Q = TRUE
)
fit2 = dsem(
  tsdata = ts(dat),
  sem = sem,
  control = control,
  family = c("fixed", "bernoulli")
)

Glm = glm(
  data = dat,
  formula = Y ~ X,
  family = binomial()
)

subset(summary(fit), path == "X -> Y")
subset(summary(fit0), path == "X -> Y")
subset(summary(fit2), path == "X -> Y")
summary(Glm)$coef['X',]

################
# 
# USE-CASE
#  Determinstic composite variable with family = "fixed" and missing values
#
################

L = lower.tri(diag(2), diag = TRUE)
X = mvtnorm::rmvnorm( n = 100, sigma = L %*% t(L) )

dat = data.frame(
  X = make_ar( X[,1], rho = 0.8),
  Y = make_ar( X[,2], rho = 0.8),
  Z = NA
)
dat[seq(1,100,by=2),'X'] = NA
dat[seq(2,100,by=2),'Y'] = NA

sem = "
  X -> X, 1, rho_X
  Y -> Y, 1, rho_Y
  X -> Z, 0, NA, 0.5
  Y -> Z, 0, NA, 0.5
  Z <-> Z, 0, NA, 0
"

fit = dsem(
  tsdata = ts(dat),
  sem = sem,
  control = dsem_control(
    gmrf = "mvn_project"
  )
)
matplot(predict(fit), type="l", lty = "solid")
matplot(dat, type = "p", add = TRUE)
