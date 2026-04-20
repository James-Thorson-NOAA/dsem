

# R = (a S) / (1 + b S) * exp(e)
# S/R = (S + b S^2) / (a S) * exp(-e)
#     = ( (1/a) + (b/a) S ) * exp(-e)
# S/R ~ lognormal( mu, sigma^2 )
#  where mu = log( (1/a) + (b/a) S )


R0 = 1
S0 = 1
h = 0.8
n = 40

SR = function(S, type = "steepness"){
  if(type=="steepness"){
    4 * h * R0 * S / (S0 * (1-h) + S * (5*h-1))
  }else if(type=="steepness"){
    a * S / ( 1 + b * S )
  }
}

S_t = runif(n)
R_t = SR(S_t) * exp(0.6*rnorm(n) - 0.6^2/2)

#devtools::install_github("James-thorson-NOAA/dsem@dev_add-family")
library(dsem)

fit = dsem(
  data = ts(data.frame(SoverR = S_t/R_t, S = S_t)),
  sem = "S -> SoverR, 0, BoverA",
  family = lognormal( )
)

