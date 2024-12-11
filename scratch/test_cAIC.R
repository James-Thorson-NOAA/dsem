
if( FALSE ){
  setwd( R'(C:\Users\James.Thorson\Desktop\Git\dsem\src\)' )
  TMB::compile( 'dsem.cpp' )
  devtools::document( R'(C:\Users\James.Thorson\Desktop\Git\dsem)' )
  devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\dsem)', dep=FALSE, force=TRUE )
}

#################
# DFA example
#################

library(dsem)
library(MARSS)
library(ggplot2)
data( harborSealWA, package="MARSS")

# Define helper function
grab = \(x,name) x[which(names(x)==name)]

# Define number of factors
# n_factors >= 3 doesn't seem to converge using DSEM or MARSS without penalties
n_factors = 2

# Add factors to data
tsdata = harborSealWA[,c("SJI","EBays","SJF","PSnd","HC")]
newcols = array( NA,
                 dim = c(nrow(tsdata),n_factors),
                 dimnames = list(NULL,paste0("F",seq_len(n_factors))) )
tsdata = ts( cbind(tsdata, newcols), start=1978)

# Scale and center (matches with MARSS does when fitting a DFA)
tsdata = scale( tsdata, center=TRUE, scale=TRUE )

#
sem = make_dfa( variables = c("SJI","EBays","SJF","PSnd","HC"),
                n_factors = n_factors )

# Initial fit
mydsem0 = dsem( tsdata = tsdata,
               sem = sem,
               family = c( rep("normal",5), rep("fixed",n_factors) ),
               estimate_delta0 = TRUE,
               control = dsem_control( quiet = TRUE,
                                       run_model = FALSE,
                                       gmrf_parameterization = "projection" ) )

# fix all measurement errors at diagonal and equal
map = mydsem0$tmb_inputs$map
map$lnsigma_j = factor( rep(1,ncol(tsdata)) )

# Fix factors to have initial value, and variables to not
map$delta0_j = factor( c(rep(NA,ncol(harborSealWA)-1), 1:n_factors) )

# Fix variables to have no stationary mean except what's predicted by initial value
map$mu_j = factor( rep(NA,ncol(tsdata)) )

# profile "delta0_j" to match MARSS (which treats initial condition as unpenalized random effect)
mydfa = dsem( tsdata = tsdata,
               sem = sem,
               family = c( rep("normal",5), rep("fixed",n_factors) ),
               estimate_delta0 = TRUE,
               control = dsem_control( quiet = TRUE,
                                       map = map,
                                       use_REML = TRUE,
                                       #profile = "delta0_j",
                                       gmrf_parameterization = "projection" ) )

cAIC(mydfa)
cAIC(mydfa, what="EDF")

###############
# Klein example
###############

# Define model
sem = "
  # Link, lag, param_name
  cprofits -> consumption, 0, a1
  cprofits -> consumption, 1, a2
  pwage -> consumption, 0, a3
  gwage -> consumption, 0, a3
  cprofits -> invest, 0, b1
  cprofits -> invest, 1, b2
  capital -> invest, 0, b3
  gnp -> pwage, 0, c2
  gnp -> pwage, 1, c3
  time -> pwage, 0, c1
"

# Load data
data(KleinI, package="AER")
TS = ts(data.frame(KleinI, "time"=time(KleinI) - 1931))
tsdata = TS[,c("time","gnp","pwage","cprofits",'consumption',
               "gwage","invest","capital")]

#
n_missing = 20
df_missing = expand.grid( seq_len(nrow(tsdata)), seq_len(ncol(tsdata)) )
df_missing = df_missing[ sample(seq_len(nrow(df_missing)), size=n_missing, replace=FALSE), ]
tsdata[ as.matrix(df_missing) ] = NA

# Fit model
fit = dsem( sem=sem,
            tsdata = tsdata,
            estimate_delta0 = TRUE,
            control = dsem_control(quiet=TRUE,
                                   getsd = FALSE,
                                   extra_convergence_checks = FALSE,
                                   newton_loops = 0) )

cAIC(fit)
cAIC(fit, what="EDF")

###################
# Linear model
###################

# simulate normal distribution
x = rnorm(100)
y = 1 + 0.5 * x + rnorm(100)
data = data.frame(x=x, y=y)

sem = "
  x -> y, 0, beta
"

# Fit as DSEM
fit = dsem( sem = sem,
            tsdata = ts(data),
            #family = c("fixed","normal"),
            control = dsem_control(quiet=TRUE) ) # gmrf_parameterization = "projection",
