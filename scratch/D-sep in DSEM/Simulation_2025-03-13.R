
###################
# ADD AXIS FOR OPTION FOR NON-NORMAL DATA
###################

# devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\dsem)', force=TRUE, dep = FALSE )

library(dsem)

root_dir = R'(C:\Users\James.Thorson\Desktop\Work files\Collaborations\2025 -- d-sep in DSEM)'
source( file.path(root_dir,"plot_histogram.R") )

# by_test and single seem to result in similar results
impute_data = c("none", "by_test", "single")[2]

# do_run_without_imputing = TRUE is very slow, and not used in paper
do_run_without_imputing = FALSE

#Date = Sys.Date()
Date = "2025-03-07"
  date_dir = file.path(root_dir, paste0(Date,"_impute=",impute_data))
  dir.create(date_dir)

if( "settings.RDS" %in% list.files(date_dir) ){
  settings = readRDS( file = file.path(date_dir,"settings.RDS") )
  attach(settings)
}else{
  # Models to run
  model_list = list(
    sem = c(
      true = "
        a -> b, 0, beta_ab, 0.5
        a -> c, 0, beta_ac, 0.5
        b -> d, 0, beta_bd, -0.5
        c -> d, 0, beta_cd, 1.0
      ",
      #backwards = "
      #  b -> a, 0, beta_ba
      #  c -> a, 0, beta_ca
      #  d -> b, 0, beta_db
      #  d -> c, 0, beta_dc
      #",
      flat = "
        a -> b, 0, beta_ab
        a -> c, 0, beta_ac
        a -> d, 0, beta_ad
      "
    ),

    dsem_simple = c(
      true = "
        a -> b, 0, rho_ab, 0.4
        a -> a, 1, rho_aa, 0.8
        b -> b, 1, rho_bb, 0.4
      ",
      #backwards = "
      #  b -> a, 0, rho_ba
      #  a -> a, 1, rho_aa
      #  b -> b, 1, rho_bb
      #",
      lagged = "
        a -> b, 1, beta_ab
        a -> a, 1, rho_aa
        b -> b, 1, rho_bb
      "
    ),

    dsem_complex = c(
      true = "
        a -> b, 0, beta_ab, 0.4
        a -> c, 0, beta_ac, 0.4
        b -> d, 0, beta_bd, 0.4
        c -> d, 0, beta_cd, 0.4
        a -> a, 1, rho_aa, 0.8
        b -> b, 1, rho_bb, 0.4
        c -> c, 1, rho_cc, 0.8
        d -> d, 1, rho_dd, 0.4
      ",
      flat = "
        a -> b, 1, beta_ab, 0.4
        a -> c, 1, beta_ac, 0.4
        a -> d, 1, beta_ad, 0.4
        a -> a, 1, rho_aa, 0.8
        b -> b, 1, rho_bb, 0.4
        c -> c, 1, rho_cc, 0.8
        d -> d, 1, rho_dd, 0.4
      "
      #backwards = "
      #  b -> a, 0, beta_ab, 0.4
      #  c -> a, 0, beta_ac, 0.4
      #  d -> b, 0, beta_bd, 0.4
      #  d -> c, 0, beta_cd, 0.4
      #  a -> a, 1, rho_aa, 0.8
      #  b -> b, 1, rho_bb, 0.4
      #  c -> c, 1, rho_cc, 0.8
      #  d -> d, 1, rho_dd, 0.4
      #"
    )
  )
  vars_set = list( letters[1:4], letters[1:2], letters[1:4] )

  #
  ntime_set = c( 25, 50, 100 )
  pmissing_set = c(0, 0.1, 0.2, 0.35, 0.5)
  n_replicates = 500

  settings = list(
    model_list = model_list,
    vars_set = vars_set,
    ntime_set = ntime_set,
    pmissing_set = pmissing_set,
    n_replicates = n_replicates
  )
  saveRDS( settings, file = file.path(date_dir,"settings.RDS") )
  capture.output(settings, file = file.path(date_dir,"settings.txt") )
}

n_scenario = length( model_list )
n_EM = max( sapply(model_list,length) )

if( "p_seplrz.RDS" %in% list.files(date_dir) ){
  p_seplrz = readRDS( file = file.path(date_dir,"p_seplrz.RDS") )
}else{
  p_seplrz = array( NA,
                  dim = c(n_scenario, n_EM, length(pmissing_set), length(ntime_set), n_replicates, 2),
                  dimnames = list( model = names(model_list),
                                   EM = 1:n_EM,
                                   pmissing = pmissing_set,
                                   n_time = ntime_set,
                                   rep = 1:n_replicates,
                                   impute = c("F","T")) )

  s_i = p_i = l_i = r_i = e_i = 1
  for( r_i in seq_len(n_replicates) ){
  for( s_i in seq_len(n_scenario) ){
  for( p_i in seq_along(pmissing_set) ){
  for( l_i in seq_along(ntime_set) ){
    if( any(is.na(p_seplrz[s_i,,p_i,l_i,r_i,])) ){
      #message("Running ", s_i, " ", o_i, " ", p_i, " ", r_i )
      vars = vars_set[[s_i]]
      p_missing = pmissing_set[p_i]
      models = model_list[[s_i]]
      n_time = ntime_set[l_i]
      message("model=", names(model_list)[s_i], " p_missing=", p_missing, " n_time=", n_time," rep=", r_i )

      # create data dims and missingness
      set.seed( r_i ) # Set seed for missingness
      data = array( 0, dim=c(n_time,length(vars)), dimnames=list(NULL,vars) )
      missing = array( rbinom(p_missing, n=prod(dim(data)), size=1), dim=dim(data), dimnames=list(NULL,vars) )
      data = ifelse(  missing == 1, NA, data )

      # simulate data
      fit0 = dsem(
        tsdata = ts(data),
        sem = models['true'],
        control = dsem_control( run_model = FALSE,
                                quiet = TRUE )
      )
      data = simulate( fit0,
                       resimulate_gmrf = TRUE,
                       seed = r_i )[[1]]  # Set seed for GMRF

      for( e_i in seq_len(n_EM) ){
        # Initial fit
        fit1 = dsem(
          tsdata = ts(data),
          sem = models[e_i],
          control = dsem_control( use_REML = FALSE,
                                  quiet = TRUE,
                                  getsd = FALSE,
                                  newton_loops = 0,
                                  extra_convergence_checks = FALSE ),
          estimate_delta0 = FALSE
        )

        # Run d-sep test
        if( fit1$opt$message == "relative convergence (4)" ){
          # impute = FALSE
          if( isTRUE(do_run_without_imputing) ){
            p_seplrz[s_i,e_i,p_i,l_i,r_i,1] = test_dsep( fit1, impute_data = "none" )
          }
          # impute = TRUE
          if(p_missing > 0){
            p_seplrz[s_i,e_i,p_i,l_i,r_i,2] = test_dsep( fit1, impute_data = impute_data )
          }else{
            p_seplrz[s_i,e_i,p_i,l_i,r_i,2] = p_seplrz[s_i,e_i,p_i,l_i,r_i,1]
          }
        }
      }
    }
  }}}}
  saveRDS( p_seplrz, file = file.path(date_dir,"p_seplrz.RDS") )
}

#################
# Plots
#################

for( p_i in seq_along(pmissing_set) ){
for( fig in c("T","F") ){
  png( file=file.path(date_dir,paste0("sim_impute=",fig,"_pmissing=",pmissing_set[p_i],".png")), width=6, height=6, res=200, units='in' )
    par( mfrow=c(n_scenario,length(ntime_set)), oma=c(4,4,4,4), mar=c(0,0,0,0), xaxs="i", yaxs="i" )
    for( s_i in seq_len(n_scenario) ){
    for( l_i in seq_along(ntime_set) ){
      p_seprz = p_seplrz[,,,l_i,,]
      by = 0.1
      plot_histogram( cbind( p_seprz[s_i,1,p_i,,fig], p_seprz[s_i,2,p_i,,"T"] ),
            breaks=seq(0,1,by=by), freq = TRUE, ylim = c(0,dim(p_seplrz)[5]) * 1.1, # ylim = c(0,1) / by * 1.1,
            xlab="", ylab = "", main = "", xaxt="n", yaxt="n",
            col = c( rgb(1,0,0,0.2), rgb(0,0,1,0.2) ) )
      if(s_i==1){
        mtext( text = ntime_set[l_i], side = 3, line=1 )
      }
      if(l_i==length(ntime_set)){
        mtext( text = names(model_list)[s_i], side = 4, line=1 )
      }
      if(s_i==n_scenario){
        axis(1)
      }
      if(l_i==1){
        axis(2)
      }
      box()
      if( s_i==1 & l_i==1 ){
        legend( "topright", fill=c( rgb(1,0,0,0.2), rgb(0,0,1,0.2) ),
                legend=c("right","wrong"), bty = "n", title="model" )
      }
      abline( h = dim(p_seplrz)[5] * by, lwd=2, lty = "dotted")
    }}
    mtext( side = 4, outer=TRUE, text = "data generating model", line=3 )
    mtext( side = 3, outer=TRUE, text = "time series length", line=2.5 )
    mtext( side = 1, outer=TRUE, text = "p-value", line=3 )
    mtext( side = 2, outer=TRUE, text = "frequency of simulations", line=2.5 )
  dev.off()
}}


for( l_i in seq_along(ntime_set) ){
for( fig in c("T","F") ){
  p_seprz = p_seplrz[,,,l_i,,]
  png( file=file.path(date_dir,paste0("sim_impute=",fig,"_length=",ntime_set[l_i],".png")), width=6, height=6, res=200, units='in' )
    par( mfrow=c(n_scenario,length(pmissing_set)), oma=c(4,4,4,4), mar=c(0,0,0,0), xaxs="i", yaxs="i" )
    for( s_i in seq_len(n_scenario) ){
    for( p_i in seq_along(pmissing_set) ){
      by = 0.1
      plot_histogram( cbind( p_seprz[s_i,1,p_i,,fig], p_seprz[s_i,2,p_i,,"T"] ),
            breaks=seq(0,1,by=by), freq = TRUE, ylim = c(0,dim(p_seplrz)[5]) * 1.1, #, ylim = c(0,1) / by * 1.1,
            xlab="", ylab = "", main = "", xaxt="n", yaxt="n",
            col = c( rgb(1,0,0,0.2), rgb(0,0,1,0.2) ) )
      if(s_i==1){
        mtext( text = pmissing_set[p_i], side = 3, line=1 )
      }
      if(p_i==length(pmissing_set)){
        mtext( text = names(model_list)[s_i], side = 4, line=1 )
      }
      if(s_i==n_scenario){
        axis(1)
      }
      if(p_i==1){
        axis(2)
      }
      box()
      if( s_i==1 & p_i==1 ){
        legend( "topright", fill=c( rgb(1,0,0,0.2), rgb(0,0,1,0.2) ),
                legend=c("right","wrong"), bty = "n", title="model" )
      }
      abline( h = dim(p_seplrz)[5] * by, lwd=2, lty = "dotted")
    }}
    mtext( side = 4, outer=TRUE, text = "data generating model", line=3 )
    mtext( side = 3, outer=TRUE, text = expression(p_missing), line=2.5 )
    mtext( side = 1, outer=TRUE, text = "p-value", line=3 )
    mtext( side = 2, outer=TRUE, text = "frequency of simulations", line=2.5 )
  dev.off()
}}


prop_seplz = apply( p_seplrz,
                    MARGIN = c(1,2,3,4,6),
                    FUN = \(v) mean(v<0.1,na.rm=TRUE) )

df = expand.grid(dimnames(prop_seplz))
df$prop = as.vector(prop_seplz)
df$EM = c("Right", "Wrong")[df$EM]
df_impute = subset( df, impute=="T" )
#df$pmissing = as.numeric(as.character(df$pmissing))

library(ggplot2)
ggplot(data=df_impute, aes(x=pmissing, y=prop, col=EM)) +
  geom_point( ) +
  facet_grid( rows = vars(model), cols = vars(n_time) ) +
  geom_hline(yintercept = 0.1) +
  labs(y= "proportion rejected", x="proportion missing data", col="Causal model") +
  ylim(0, 1) +
  #theme_classic( )
  theme(panel.border = element_rect(color = "black", fill = NA),
       panel.background = element_blank(),
       panel.grid = element_blank(),
       panel.spacing.x = unit(0,"line"))
ggsave( file=file.path(date_dir,"Fig_4_power_analysis.png"), width=6, height=6 )

df_noimpute = subset( df, impute=="F" )
library(ggplot2)
ggplot(data=df_noimpute, aes(x=pmissing, y=prop, col=EM)) +
  geom_point( ) +
  facet_grid( rows = vars(model), cols = vars(n_time) ) +
  geom_hline(yintercept = 0.1) +
  labs(y= "proportion rejected", x="proportion missing data", col="Causal model") +
  ylim(0, 1) +
  #theme_classic( )
  theme(panel.border = element_rect(color = "black", fill = NA),
       panel.background = element_blank(),
       panel.grid = element_blank(),
       panel.spacing.x = unit(0,"line"))
ggsave( file=file.path(date_dir,"Fig_4_power_analysis-noimpute.png"), width=6, height=6 )

library(ggplot2)
ggplot(data=df, aes(x=pmissing, y=prop, col=EM, size=impute) ) +
  geom_point( ) +
  facet_grid( rows = vars(model), cols = vars(n_time) ) +
  geom_hline(yintercept = 0.1) +
  labs(y= "proportion rejected", x="proportion missing data", col="Causal model") +
  ylim(0, 1) +
  #theme_classic( )
  theme(panel.border = element_rect(color = "black", fill = NA),
       panel.background = element_blank(),
       panel.grid = element_blank(),
       panel.spacing.x = unit(0,"line"))
ggsave( file=file.path(date_dir,"Fig_4_power_analysis-both.png"), width=6, height=6 )

##############
# Show graphs
##############

library(igraph)

s_i = e_i = 1
png( file=file.path(date_dir,paste0("Fig_1_graphs.png")), width=4, height=6, res=200, units='in' )
  par( mfrow=c(n_scenario,n_EM), mar=c(0,0,0,0), oma=c(0,2,2,0) )
  for( s_i in seq_len(n_scenario) ){
  for( e_i in seq_len(n_EM) ){
    vars = vars_set[[s_i]]
    models = model_list[[s_i]]
    n_time = ntime_set[1]
    data = array( 0, dim=c(n_time,length(vars)), dimnames=list(NULL,vars) )
    if( s_i==1 ){
      layout = rbind( a = c(0,1), b = c(-1,0), c=c(1,0), d=c(0,-1) )
    }
    if( s_i==2 ){
      layout = rbind( a = c(-1,1), b = c(1,-1) )
    }
    if( s_i==3 ){
      layout = rbind( a = c(0,1), b = c(-1,0), c=c(1,0), d=c(0,-1) )
    }

    # simulate data
    fit0 = dsem(
      tsdata = ts(data),
      sem = models[e_i],
      control = dsem_control( run_model = FALSE,
                              quiet = TRUE )
    )
    df = fit0$sem_full
    df = subset( df, direction==1 )
    df_graph = data.frame(from = df$first, to = df$second, label = "" ) # df$start
    pg <- graph_from_data_frame(d = df_graph, directed=TRUE, vertices=vars )
    #pg = simplify(pg)
    plot( pg, vertex.shape="rectangle", vertex.size=0, vertex.size2=0, vertex.label.cex=2,
          vertex.color="grey", vertex.label.color="black", edge.label.color="black",
          edge.label.cex=1.5, xlim=1.6*c(-1,1), ylim=1.6*c(-1,1),
          edge.color = ifelse(df$lag==1, "blue", "black"), loop.size = 2,
          arrow.size = 0.8, layout = layout )  # layout=layout,
    if(s_i==1){
      mtext( text = c("Right","Wrong")[e_i], side = 3, line=0 )
    }
    if(e_i==1){
      mtext( text = names(model_list)[s_i], side = 2, line=0 )
    }
  }}
dev.off()

# Conditional independence for DSEM_simple
s_i = 2
e_i = 1
vars = vars_set[[s_i]]
models = model_list[[s_i]]
n_time = ntime_set[1]
data = array( 0, dim=c(n_time,length(vars)), dimnames=list(NULL,vars) )

# simulate data
fit0 = dsem(
  tsdata = ts(data),
  sem = models[e_i],
  control = dsem_control( run_model = FALSE,
                          quiet = TRUE )
)
data = simulate( fit0,
                 resimulate_gmrf = TRUE,
                 seed = 1 )[[1]]  # Set seed for GMRF
fit1 = dsem(
  tsdata = ts(data),
  sem = models[e_i],
  control = dsem_control( use_REML = FALSE,
                          quiet = TRUE,
                          getsd = FALSE,
                          newton_loops = 0,
                          extra_convergence_checks = FALSE ),
  estimate_delta0 = FALSE
)
order = expand.grid(vars, paste0("lag",0:2))
order = apply( order, MARGIN=1, FUN=paste0, collapse="_")
out = test_dsep(fit1, what="all", order = order)

png( file=file.path(date_dir,paste0("Fig_2_conditional_independence_tests.png")), width=3, height=4, res=200, units='in' )
  nrows = ceiling( sqrt(length(out$sems)))
  par( mfrow=c(nrows,ceiling(length(out$sems)/nrows)), mar=c(0,0,2,0), oma=c(0,0,0,0) )
  for( z_i in seq_along(out$sems) ){
    df = summary(out$fits[[z_i]])
    df = subset( df, direction==1 )
    df_graph = data.frame(from = df$first, to = df$second, label = "" ) # df$start
    pg <- graph_from_data_frame(d = df_graph, directed=TRUE, vertices=vars )
    #pg = simplify(pg)
    plot( pg, vertex.shape="rectangle", vertex.size=0, vertex.size2=0, vertex.label.cex=2,
          vertex.color="grey", vertex.label.color="black", edge.label.color="black",
          edge.label.cex=1.5, xlim=1.6*c(-1,1), ylim=1.6*c(-1,1),
          edge.color = c("black", "blue", "green")[df$lag+1], loop.size = 2,
          arrow.size = 0.6, edge.lty=c("solid",rep("dotted",nrow(df)-1)) )  # layout=layout,
    title( paste0("Conditional independence ", z_i) )
    if(z_i == 1){
      legend( "topleft", bty="n", fill = c("black", "blue", "green"),
            legend = 0:2, title = "lag" )
    }
  }
dev.off()

###################
# Verbal example
###################

vars = letters[1:3]
A = array(0, dim=rep(length(vars),2), dimnames=list(vars,vars))
A[cbind(c("a","b"),c("b","c"))] = 1
ggm::basiSet(A)

