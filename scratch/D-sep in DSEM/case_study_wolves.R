###############
# Wolf-Moose
###############


setwd( R'(C:\Users\James.Thorson\Desktop\Work files\Collaborations\2025 -- d-sep in DSEM)')

library(dsem)

data(isle_royale)
data = ts( log(isle_royale[,2:3]), start=1959)
colnames(data) = c("W", "M")

# Density dependence
sem1 = c(
  "W -> W, 1, rho1_ww",
  "M -> M, 1, rho1_mm"
)
fit1 = dsem(
  tsdata = data,
  sem = sem1,
  estimate_delta0 = TRUE
)
AIC1 = AIC(fit1)

# Bottom-up
sem2 = c(
  "W -> W, 1, rho1_ww",
  "M -> W, 1, rho1_mw",
  "M -> M, 1, rho1_mm"
)
fit2 = dsem(
  tsdata = data,
  sem = sem2,
  estimate_delta0 = TRUE
)
AIC2 = AIC(fit2)

# Top-down
sem3 = c(
  # Lag-1
  "W -> W, 1, rho1_ww",
  "W -> M, 1, rho1_wm",
  "M -> M, 1, rho1_mm"
)
fit3 = dsem(
  tsdata = data,
  sem = sem3,
  estimate_delta0 = TRUE
)
AIC3 = AIC(fit3)

# Both
sem4 = c(
  # Lag-1
  "W -> W, 1, rho1_ww",
  "M -> W, 1, rho1_mw",
  "W -> M, 1, rho1_wm",
  "M -> M, 1, rho1_mm"
)
fit4 = dsem(
  tsdata = data,
  sem = sem4,
  estimate_delta0 = TRUE
)
AIC4 = AIC(fit4)

library(igraph)
layout = cbind( x=c(1, 1, 2, 2), y=c(1,2,1,2) )
  var_labels = c( "M(t)", "W(t)", "M(t+1)", "W(t+1)" )
  rownames(layout) = var_labels
  config_labels = c("Density dependence", "Bottom up", "Top down", "Both")
png( file=file.path(getwd(),"Fig_6_graphs.png"), width=6, height=6, res=200, units="in" )
  par( mfrow=c(2,2), mar=c(1,0,2,0) )
  for( i in seq_along(config_labels) ){
    fit = list( fit1, fit2, fit3, fit4 )[[i]]
    test = test_dsep(fit)
    AIC_i = c( AIC1, AIC2, AIC3, AIC4 )[i] - min(c(AIC1,AIC2,AIC3,AIC4))
    df = summary(fit)
    df = subset( df, direction==1 )
    df = data.frame(from = paste0(df$first,"(t)"), to = paste0(df$second,"(t+1)"), label = round(df$Estimate,2) )
    pg <- graph_from_data_frame(d = df, directed=TRUE, vertices=data.frame(var_labels) )
    #pg = simplify(pg)
    plot( pg, vertex.shape="rectangle", vertex.size=0, vertex.size2=0, vertex.label.cex=2,
          vertex.color="grey", vertex.label.color="black", edge.label.color="black",
          edge.label.cex=1.5, layout=layout, xlim=1.2*c(-1,1), ylim=1.2*c(-1,1) + c(0,0.5) )
    title( config_labels[i] )
    #mtext( side=1, text=round(AIC_table[i,2],1) )
    legend( "topright", bty="n", legend=paste0("?AIC = ",round(AIC_i,1)) )
    legend( "topleft", bty="n", legend=paste0("p = ", prettyNum(test, digits=2)) )
    box()
  }
dev.off()


## AIC selection
#selexAIC = stepwise_selection(
#  model_options = sem,
#  model_shared = "",
#  tsdata = data,
#  estimate_delta0 = TRUE
#)
#fitAIC = dsem( sem = selexAIC$model,
#           tsdata = data,
#           estimate_delta0 = TRUE )
#summary(fitAIC)
#test_dsep( fitAIC )

