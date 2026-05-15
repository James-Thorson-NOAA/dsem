
library(dsem)

# Settings
n_factors = 5
n_vars = 10
n_t = 50

# Simulate factors as random-walk
z_tf = matrix(rnorm(n_t*n_factors), ncol = n_factors)
z_tf = apply(z_tf, MARGIN = 2, FUN = cumsum)

# Could rotate or orthogonalize simulated factors
# BUT not doing it to keep things simple

# Simulate loadings
L_cf = matrix(rnorm(n_vars*n_factors), ncol = n_factors)
# lower triangle here to simplify comparison of estimated and true
L_cf[upper.tri(L_cf)] = 0

#
x_tc = z_tf %*% t(L_cf)

# Measurements
y_tc = x_tc + 0.1 * matrix(rnorm(n_vars*n_t),ncol=n_vars)
colnames(y_tc) = letters[seq_len(n_vars)]

# Factor block
F_tf = matrix( NA, nrow = n_t, ncol = n_factors )
colnames(F_tf) = paste0("F",seq_len(n_factors))

#
sem = make_dfa(
  variables = colnames(y_tc),
  n_factors = n_factors,
  factor_names = colnames(F_tf)
)

#
fit = dsem(
  sem = sem,
  family = c(rep("normal", n_vars),rep("fixed",n_factors)),
  tsdata = ts(cbind(y_tc,F_tf)),
  control = dsem_control(
    gmrf_parameterization = "projection",
    trace = 1
  )
)

# Extract loadings
Lhat_cf = matrix(0, ncol = n_factors, nrow = n_vars)
Lhat_cf[lower.tri(Lhat_cf,diag=TRUE)] = subset(summary(fit),direction==1 & lag==0)$Estimate

# Deal with label-swiching signs
Lhat_cf = sweep(
  Lhat_cf,
  MARGIN = 2,
  STAT = sign(diag(cov(L_cf, Lhat_cf))),
  FUN = "*"
)

# Compare
plot(
  x = L_cf,
  y = Lhat_cf
)
