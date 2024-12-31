
library(dsem)
library(RTMB)

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

# Fit model
fit = dsem( sem=sem,
            tsdata = tsdata,
            estimate_delta0 = TRUE,
            control = dsem_control( quiet=TRUE ) )
#ParHat = fit$internal$parhat
#Rep = fit$obj$report()

# RUN dsemRTMB line-by-line
if( FALSE ){
  control = dsem_control()
  covs = colnames(tsdata)
}

#
source( file.path(R'(C:\Users\James.Thorson\Desktop\Git\dsem\R)', "make_matrices.R") )
source( file.path(R'(C:\Users\James.Thorson\Desktop\Git\dsem\R)', "get_jnll.R") )
source( file.path(R'(C:\Users\James.Thorson\Desktop\Git\dsem\R)', "read_model.R") )
source( file.path(R'(C:\Users\James.Thorson\Desktop\Git\dsem\R)', "dsemRTMB.R") )
fitRTMB = dsemRTMB( sem = sem,
            tsdata = tsdata,
            estimate_delta0 = TRUE,
            control = dsem_control( quiet = TRUE) )

