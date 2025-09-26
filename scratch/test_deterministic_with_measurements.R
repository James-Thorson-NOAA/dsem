
# devtools::install_local( R'(C:\Users\james\OneDrive\Desktop\Git\dsem)', force = TRUE, dep = FALSE )

library(dsem)

################
# 
# USE-CASE
#  Determinstic and family = "fixed" with measurements
#
################


p_missing = 0.2

X = cumsum( rnorm(100) )
Y = 2 * X
Z = 0.5*Y + rnorm(100)
X_missing = sample(seq_along(X), size=p_missing*length(X), replace=FALSE)

# Missing-ness
Xobs = X
Xobs[X_missing] = NA
Yobs = Y
Yobs[-X_missing] = NA

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
for(t in 2:length(X)) X[t] = 0.8 * X[t-1] + sqrt(1-0.8) * X[t]
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

control = dsem_control(
  gmrf_parameterization = "mvn_project",
  use_REML = FALSE
)
fit = dsem(
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
summary(Glm)$coef['X',]
