

library(dsem)

##########
# Linear model example
##########

x = rnorm(100)
y = 1 + 1 * x + rnorm(100)
data = data.frame(x=x, y=y)

# Fit as linear model
Lm = lm( y ~ x, data=data )

# Fit as DSEM
fit = dsem( sem = "x -> y, 0, beta",
            tsdata = ts(data),
            control = dsem_control(quiet=TRUE) )

#
partition_variance( fit,
                    which_pred = "x",
                    which_response = "y",
                    n_lags = 10 )

#################
# Klein example
#################
data(KleinI, package="AER")
TS = ts(data.frame(KleinI, "time"=time(KleinI) - 1931))

# Specify by declaring each arrow and lag
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
tsdata = TS[,c("time","gnp","pwage","cprofits",'consumption',
               "gwage","invest","capital")]
fit = dsem( sem=sem,
            tsdata = tsdata,
            estimate_delta0 = TRUE,
            control = dsem_control(
              quiet = TRUE,
              newton_loops = 0) )

#
partition_variance( fit,
                    which_pred = "gnp",
                    which_response = "pwage",
                    n_lags = 10 )
