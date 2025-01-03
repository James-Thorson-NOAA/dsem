

library(dsem)
library(RTMB)
library(Matrix)

# Files
source( file.path(R'(C:\Users\James.Thorson\Desktop\Git\dsem\R)', "make_matrices.R") )
source( file.path(R'(C:\Users\James.Thorson\Desktop\Git\dsem\R)', "compute_nll.R") )
source( file.path(R'(C:\Users\James.Thorson\Desktop\Git\dsem\R)', "read_model.R") )
source( file.path(R'(C:\Users\James.Thorson\Desktop\Git\dsem\R)', "dsemRTMB.R") )

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
fit = dsemRTMB( sem=sem,
            tsdata = tsdata,
            estimate_delta0 = TRUE,
            #family = rep("normal",ncol(tsdata)),
            control = dsem_control( quiet = FALSE,
                                    #run_model = TRUE,
                                    #use_REML = TRUE,
                                    gmrf_parameterization = "separable" ) )

fit$obj$simulate()$y_tj
