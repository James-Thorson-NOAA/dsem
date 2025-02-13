
# devtools::document( R'(C:\Users\James.Thorson\Desktop\Git\dsem)' )
# devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\dsem)', force=TRUE, dep=FALSE )
# devtools::install_github( "James-thorson-NOAA/dsem@dev", force=TRUE, dep=FALSE )

library(dsem)
library(ggm)
library(igraph)


##############
# No time simulation
##############

# Confirm that graph has conditional dependencies with non-null conditioning set
A <- DAG(x5~ x3+x4, x3~ x2, x4~x2)
basiSet(A)
plot(graph_from_adjacency_matrix(A))

# Settings
p_missing = 0.25
use_imputed_data = FALSE

# Simulation loop
pvalue_rz = array(NA, dim=c(100,4) )
for( r in seq_len(dim(pvalue_rz)[1]) ){
  if( any(is.na(pvalue_rz[r,])) ){
    set.seed(r)

    # simulate normal distribution
    a = rnorm(100)
    b = 1 + 0.5 * a + rnorm(100)
    c = 1 + 0.5 * a + rnorm(100)
    d = 2 + -0.5*b + 1*c + rnorm(100)
    data = data.frame(a=a, b=b, c=c, d=d)

    # Missing at random
    missing = array( rbinom(p_missing, n=prod(dim(data)), size=1), dim=dim(data))
    data = ifelse( missing == 1, NA, as.matrix(data) )
    colnames( data ) = letters[1:4]

    # Correct SEM
    sem1 = "
      a -> b, 0, beta_ab
      a -> c, 0, beta_ac
      b -> d, 0, beta_bd
      c -> d, 0, beta_cd
    "
    fit1 = dsem(
      tsdata = ts(data),
      sem = sem1
    )
    test1 = test_dsep( fit1, use_imputed_data = FALSE, what = "all" )
    test2 = test_dsep( fit1, use_imputed_data = TRUE, what = "all" )
    pvalue_rz[r,1] = test1$pvalue
    pvalue_rz[r,2] = test2$pvalue

    # Backwards causality
    sem2 = "
      b -> a, 0, beta_ba
      c -> a, 0, beta_ca
      d -> b, 0, beta_db
      d -> c, 0, beta_dc
    "
    fit2 = dsem(
      tsdata = ts(data),
      sem = sem2
    )
    pvalue_rz[r,3] = test_dsep( fit2 )

    # Flat regression structure
    sem3 = "
      a -> b, 0, beta_ab
      a -> b, 0, beta_ab
      a -> c, 0, beta_ac
    "
    fit3 = dsem(
      tsdata = ts(data),
      sem = sem2
    )
    pvalue_rz[r,4] = test_dsep( fit3 )
  }
}

# Show tests
par( mfrow=c(2,2), oma=c(2,0,0,0) )
hist(pvalue_rz[,1], breaks=seq(0,1,by=0.05), main="right model")
hist(pvalue_rz[,2], breaks=seq(0,1,by=0.05), main="right (imputed data)")
hist(pvalue_rz[,3], breaks=seq(0,1,by=0.05), main="backwards causality")
hist(pvalue_rz[,4], breaks=seq(0,1,by=0.05), main="regression-style structure")
mtext( side=1, text="p-value (rejecting test that model is correct)", outer=TRUE )

######### stepwise selection on both directions
model_options = c(
  "a -> b, 0, beta_ab",
  "a -> c, 0, beta_ac",
  "a -> d, 0, beta_ad",
  "b -> c, 0, beta_bc",
  "b -> d, 0, beta_bd",
  "c -> d, 0, beta_cd",

  "b -> a, 0, beta_ba",
  "c -> a, 0, beta_ca",
  "d -> a, 0, beta_da",
  "c -> b, 0, beta_cb",
  "d -> b, 0, beta_db",
  "d -> c, 0, beta_dc"
)
selex = stepwise_selection(
  model_options = model_options,
  model_shared = "",
  tsdata = ts(data)
)
cat(selex$model)

fit_selex = dsem(
  sem = selex$model,
  tsdata = ts(data)
)

#################
# VAR simulation
#################

model = c( "simple", "complex" )[2]
p_missing = 0.5

if( model == "simple" ){
  vars = letters[1:2]
  sem0 = "
    a -> b, 0, rho_ab, 0.4
    a -> a, 1, rho_aa, 0.8
    b -> b, 1, rho_bb, 0.4
  "
  sem1 = sem2 = sem0
  sem3 = "
    b -> a, 0, rho_ba
    a -> a, 1, rho_aa
    b -> b, 1, rho_bb
  "
  sem4 = "
    a -> b, 1, beta_ab
    a -> a, 1, rho_aa
    b -> b, 1, rho_bb
  "
}
if( model == "complex" ){
  vars = letters[1:4]
  sem0 = "
    a -> b, 0, beta_ab, 0.4
    a -> c, 0, beta_ac, 0.4
    b -> d, 0, beta_bd, 0.4
    c -> d, 0, beta_cd, 0.4
    a -> a, 1, rho_aa, 0.8
    b -> b, 1, rho_bb, 0.4
    c -> c, 1, rho_cc, 0.8
    d -> d, 1, rho_dd, 0.4
  "
  sem1 = sem2 = sem0
  #sem3 = "
  #  a -> b, 1, beta_ab, 0.4
  #  a -> c, 1, beta_ac, 0.4
  #  b -> d, 1, beta_bd, 0.4
  #  c -> d, 1, beta_cd, 0.4
  #  a -> a, 1, rho_aa, 0.8
  #  b -> b, 1, rho_bb, 0.4
  #  c -> c, 1, rho_cc, 0.8
  #  d -> d, 1, rho_dd, 0.4
  #"
  sem3 = "
    b -> a, 0, beta_ab, 0.4
    c -> a, 0, beta_ac, 0.4
    d -> b, 0, beta_bd, 0.4
    d -> c, 0, beta_cd, 0.4
    a -> a, 1, rho_aa, 0.8
    b -> b, 1, rho_bb, 0.4
    c -> c, 1, rho_cc, 0.8
    d -> d, 1, rho_dd, 0.4
  "
}

# Simulation loop
n_replicates = 500
n_t = 100
pvalue_rz = array(NA, dim=c(n_replicates,4) )
for( r in seq_len(n_replicates) ){

  if( any(is.na(pvalue_rz[r,1:4])) ){
    set.seed(r)
    message("Running ", r)

    # create data dims and missingness
    data = array( 0, dim=c(100,length(vars)), dimnames=list(NULL,vars) )
    missing = array( rbinom(p_missing, n=prod(dim(data)), size=1), dim=dim(data), dimnames=list(NULL,vars) )
    data = ifelse(  missing == 1, NA, data )

    # simulate data
    fit0 = dsem(
      tsdata = ts(data),
      sem = sem0,
      control = dsem_control( run_model = FALSE,
                              quiet = TRUE )
    )
    data = simulate( fit0,
                     resimulate_gmrf = TRUE )[[1]]

    # Correct SEM
    fit1 = dsem(
      tsdata = ts(data),
      sem = sem1,
      control = dsem_control( use_REML = FALSE,
                              quiet = TRUE ),
      estimate_delta0 = FALSE
    )
    test1 = test_dsep( fit1,
      #test = "lr"
      what = "all" )
    # system.time( test_dsep(fit1, test="wald") )
    pvalue_rz[r,1] = test1$pvalue

    test2 = test_dsep( fit1,
      #test = "lr"
      what = "all",
      use_imputed_data = TRUE )
    pvalue_rz[r,2] = test2$pvalue

    # Wrong SEM
    fit3 = dsem(
      tsdata = ts(data),
      control = dsem_control( use_REML = FALSE,
                              quiet = TRUE ),
      sem = sem3,
      estimate_delta0 = FALSE
    )
    test3 = test_dsep( fit3,
      #test = "lr",
      what = "all" )
    pvalue_rz[r,3] = test3$pvalue

    test4 = test_dsep( fit3,
      #test = "lr"
      what = "all",
      use_imputed_data = TRUE )
    pvalue_rz[r,4] = test4$pvalue
#    # Wrong SEM
#    fit4 = dsem(
#      tsdata = ts(data),
#      control = dsem_control( use_REML = FALSE ),
#      sem = sem4,
#      estimate_delta0 = FALSE
#    )
#    test4 = test_dsep( fit4,
#      #test = "lr",
#      what = "all" )
#    pvalue_rz[r,4] = test4$pvalue
  }
}

# Show tests
par( mfrow=c(2,2), oma=c(2,0,0,0) )
hist(pvalue_rz[,1], breaks=seq(0,1,by=0.05), main="right model")
hist(pvalue_rz[,2], breaks=seq(0,1,by=0.05), main="right (imputed data)")
hist(pvalue_rz[,3], breaks=seq(0,1,by=0.05), main="backwards causality")
hist(pvalue_rz[,4], breaks=seq(0,1,by=0.05), main="backwards (imputed data)")
mtext( side=1, text="p-value (rejecting test that model is correct)", outer=TRUE )

# Check differences
sems1 = sapply( test1$sems, FUN = \(x){paste0(x,collapse=" \n ")} )
sems3 = sapply( test3$sems, FUN = \(x){paste0(x,collapse=" \n ")} )


#############
# Real data
#############

# Bering Sea
data(bering_sea)
Z = ts( bering_sea )
family = rep('fixed', ncol(bering_sea))

# Specify model
sem = "
  # Link, lag, param_name
  log_seaice -> log_CP, 0, seaice_to_CP
  log_CP -> log_Cfall, 0, CP_to_Cfall
  log_CP -> log_Esummer, 0, CP_to_E
  log_PercentEuph -> log_RperS, 0, Seuph_to_RperS
  log_PercentCop -> log_RperS, 0, Scop_to_RperS
  log_Esummer -> log_PercentEuph, 0, Esummer_to_Suph
  log_Cfall -> log_PercentCop, 0, Cfall_to_Scop
  SSB -> log_RperS, 0, SSB_to_RperS

  log_seaice -> log_seaice, 1, AR1, 0.001
  log_CP -> log_CP, 1,  AR2, 0.001
  log_Cfall -> log_Cfall, 1, AR4, 0.001
  log_Esummer -> log_Esummer, 1, AR5, 0.001
  SSB -> SSB, 1, AR6, 0.001
  log_RperS ->  log_RperS, 1, AR7, 0.001
  log_PercentEuph -> log_PercentEuph, 1, AR8, 0.001
  log_PercentCop -> log_PercentCop, 1, AR9, 0.001
"

# Fit
fit = dsem( sem = sem,
            tsdata = Z,
            family = family,
            control = dsem_control(use_REML=FALSE, quiet=TRUE) )
test = test_dsep( fit,
                  what = "all" )

# Sort by missing links
order = order( test$pvalues, decreasing = FALSE )
links = sapply( test$sem,
        FUN = function(vec) vec[1] )
data.frame( link=links, test$pvalues )[order,]

###############
# Wolf-Moose
###############

data(isle_royale)
data = ts( log(isle_royale[,2:3]), start=1959)

sem = "
  # Link, lag, param_name
  wolves -> wolves, 1, arW
  moose -> wolves, 1, MtoW
  wolves -> moose, 1, WtoM
  moose -> moose, 1, arM
  moose -> wolves, 0, cross_cor
"
# initial first without delta0 (to improve starting values)
fit = dsem( sem = sem,
             tsdata = data,
             estimate_delta0 = FALSE )
test = test_dsep( fit, what = "all" )
test$pvalue
test$CIC

