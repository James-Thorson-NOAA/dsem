
library(dsem)
#library(phylopath)

setwd( R'(C:\Users\James.Thorson\Desktop\Work files\Collaborations\2025 -- d-sep in DSEM)' )

Survey = read.csv("expected_shelikof_indices_Jan2022.csv")
Covs = read.csv("CatchabilityCovariates_2023-03-29.csv")
colnames(Covs) = c("Year", "P", "A", "T" )
Data = merge(Survey,Covs)
Data$Q = log(Data[,'Observed3p'] / Data[,'No_survey_pred'])
#
#Data$log_pred = log(Data[,'No_survey_pred'])

# Scaling to avoid numerical stiffness
Data$A = (Data$A - 30)
#Data$P = Data$P + 3
Data$T = Data$T - 3

#################
# ChECK PAPER
#  ... all look right
#################

summary(lm(Q ~ A, data=Data))
summary(lm(Q ~ P, data=Data))
summary(lm(Q ~ T, data=Data))

#################
# Without latent
#################

Data2 = cbind(Data[,c('Q','P','A','T')])     # 'log_pred',
Z = ts( Data2, start=1983 )
Z = subset(Z, time(Z) >= 1992 )
Z = ts(Z, start=1992)

sem = "
  # Links
  T -> A, 0, b1
  A -> P, 0, b3
  A -> Q, 0, b4

  # ARs
  A -> A, 1, rho1
  T -> T, 1, rho2
  P -> P, 1, rho3
  Q -> Q, 1, rho5
"

fit1 = dsem( tsdata = Z,
            sem = sem )
test1 = test_dsep(fit1)
AIC1 = AIC(fit1)

#pvalue = rep(NA, 20)
#for(i in seq_along(pvalue)) pvalue[i] = test_dsep(fit1, seed=i ) #

##########
# OTHER STRUCTURE
##########

sem = "
  # Links
  T -> A, 0, b1
  T -> P, 0, b3
  T -> Q, 0, b4

  # ARs
  A -> A, 1, rho1
  T -> T, 1, rho2
  P -> P, 1, rho3
  Q -> Q, 1, rho5
"

fit2 = dsem( tsdata = Z,
            sem = sem )
test2 = test_dsep(fit2)
AIC2 = AIC(fit2)

#pvalue = rep(NA, 20)
#for(i in seq_along(pvalue)) pvalue[i] = test_dsep(fit2, seed=i ) #

##########
# OTHER STRUCTURE
##########

sem = "
  # Links
  T -> Q, 0, b1
  A -> Q, 0, b3
  P -> Q, 0, b4

  # ARs
  A -> A, 1, rho1
  T -> T, 1, rho2
  P -> P, 1, rho3
  Q -> Q, 1, rho5
"

fit3 = dsem( tsdata = Z,
            sem = sem )
test3 = test_dsep(fit3)
AIC3 = AIC(fit3)

#pvalue = rep(NA, 20)
#for(i in seq_along(pvalue)) pvalue[i] = test_dsep(fit3, seed=i ) #

####
library(igraph)

#layout = cbind( x=c(1, 2, 2, 3), y=c(2,1,3,2) )
  #var_labels = c( "P", "Q", "T", "A" )
  #rownames(layout) = var_labels
layout = rbind( P = c(3,2), Q = c(2,1), T = c(2,3), A = c(1,2) )
  config_labels = c("Temperature as driver", "Availability regression", "Timing as mediator" )
png( file=file.path(getwd(),"Fig_5_graphs.png"), width=6, height=3, res=200, units="in" )
  par( mfrow=c(1,3), mar=c(1,0,2,0) )
  for( i in seq_along(config_labels) ){
    fit = list( fit2, fit3, fit1 )[[i]]
    test = test_dsep(fit)
    AIC_i = c( AIC1, AIC2, AIC3 )[i] - min(c(AIC1,AIC2,AIC3))
    df = summary(fit)
    df = subset( df, first!=second )
    df = data.frame(from = df$first, to = df$second, label = round(df$Estimate,2) )
    pg <- graph_from_data_frame(d = df, directed=TRUE, vertices=data.frame(rownames(layout)) )
    #pg = simplify(pg)
    plot( pg, vertex.shape="rectangle", vertex.size=0, vertex.size2=0, vertex.label.cex=2,
          vertex.color="grey", vertex.label.color="black", edge.label.color="black",
          edge.label.cex=1.5, layout=layout, xlim=1.2*c(-1,1), ylim=1.2*c(-1,1) )
    title( config_labels[i] )
    #mtext( side=1, text=round(AIC_table[i,2],1) )
    legend( "topright", bty="n", legend=paste0("?AIC = ",round(AIC_i,1)) )
    legend( "topleft", bty="n", legend=paste0("p = ", prettyNum(test, digits=2)) )
    box()
  }
dev.off()
