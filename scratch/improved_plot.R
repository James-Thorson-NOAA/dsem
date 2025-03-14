
library(dsem)
sem = "
  a -> b, 0, ab, 0.2
  b -> c, 0, bc, -0.2
  a -> a, 1, aa_lag1, 0.8
  b -> b, 1, bb_lag1, 0.6
"
dat = ts(data.frame(a=rep(0,10), b=rep(0,10), c=rep(0,10)))
out0 = dsem( tsdata = dat,
             sem = sem,
             control = dsem_control(run_model=FALSE) )
out0$internal$parhat = out0$obj$env$parList()

P_kk = dsem:::get_P( out0, times = seq_len(max(out)) )
pg = graph_from_adjacency_matrix(Matrix::t(P_kk), weighted="directed")
coords = layout_(pg, with_sugiyama())
plot( pg, layout = coords )
