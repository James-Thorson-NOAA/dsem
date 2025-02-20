

################### 
# COPIED from package ggm
# With permission from maintainer Giovanni Marchetti
# by email Jan. 19, 2025
# and under licence GPL-2
###################


basiSet <-
function(amat) {
  ### Basis set of a DAG with adjacency matrix amat.
  amat <- topSort(amat)
  nod <- rownames(amat)
  dv <- length(nod)
  ind <- NULL
  ## NOTE. This is correct if the adj mat is upper triangular.
  for (r in 1:dv) {
    for (s in r:dv) {
      if ((amat[r, s] != 0) | (s == r)) {
        next
      } else {
        ed <- nod[c(r, s)]
        pa.r <- nod[amat[, r] == 1]
        pa.s <- nod[amat[, s] == 1]
        dsep <- union(pa.r, pa.s)
        dsep <- setdiff(dsep, ed)
        b <- list(c(ed, dsep))
        ind <- c(ind, b)
      }
    }
  }
  ##      ind <- lapply(ind, function(x) nn[x])
  ind
}

findPath <-
function(amat, st, en, path = c()) {
  ### Find a path between nodes st and en in a UG with adjacency mat. amat.
  indices <- 1:nrow(amat)
  if (st == en) { # st is 'node' in recursive calls
    return(c(path, st))
  }
  if (sum(amat[st, ]) == 0) {
    return(NULL)
  }
  ## ne <- bd(st,amat)
  ne <- indices[amat[st, ] == 1] # Boundary of x. Assumes that amat is symmetric
  for (node in ne) {
    if (!is.element(node, c(path, st))) {
      newpath <- findPath(amat, node, en, c(path, st))
      if (!is.null(newpath)) {
        return(newpath)
      }
    }
  }
}

isAcyclic <-
function(amat, method = 2) {
  ### Tests if the graph is acyclic.
  if (method == 1) {
    G <- graph.adjacency(amat)
    return(max(clusters(G, mode = "strong")$csize) == 1)
  } else if (method == 2) {
    B <- transClos(amat)
    l <- B[lower.tri(B)]
    u <- t(B)[lower.tri(t(B))]
    com <- (l & u)
    return(all(!com))
  } else {
    stop("Wrong method.")
  }
}

transClos <-
function(amat) {
  ### Transitive closure of the relation with adjacency matrix amat.
  if (nrow(amat) == 1) {
    return(amat)
  }
  A <- amat
  diag(A) <- 1
  repeat {
    B <- sign(A %*% A)
    if (all(B == A)) {
      break
    } else {
      A <- B
    }
  }
  diag(A) <- 0
  A
}

topOrder <-
function(amat) {
  ### Return the nodes in topological order (parents before children).
  ### Translated from: Kevin Murphy's BNT.
  if (!isAcyclic(amat)) stop("The graph is not acyclic!")
  n <- nrow(amat)
  nod <- 1:n
  indeg <- rep(0, n)
  up <- !amat[lower.tri(amat)]
  if (all(up)) {
    return(nod)
  }
  zero.indeg <- c() #  a stack of nodes with no parents
  for (i in nod) {
    indeg[i] <- sum(amat[, i])
    if (indeg[i] == 0) {
      zero.indeg <- c(i, zero.indeg)
    }
  }
  s <- 1
  ord <- rep(0, n)
  while (length(zero.indeg) > 0) {
    v <- zero.indeg[1] #  pop v
    zero.indeg <- zero.indeg[-1]
    ord[s] <- v
    s <- s + 1
    cs <- nod[amat[v, ] == 1]
    if (length(cs) == 0) next
    for (j in 1:length(cs)) {
      k <- cs[j]
      indeg[k] <- indeg[k] - 1
      if (indeg[k] == 0) {
        zero.indeg <- c(k, zero.indeg)
      } # push k
    }
  }
  ord
}

topSort <-
function(amat) {
  ### Topological sort of the DAG with adjacency matrix amat.
  ord <- topOrder(amat)
  amat[ord, ord]
}
