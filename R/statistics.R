#' Compute the probability of vertex i connecting to vertex j.
#'
#' @param ki integer degree of vertex i
#' @param kj integer degree of vertex j
#' @param M integer size of graph G, the total number of edges within G.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

comp_pij <- function(ki, kj, M, sparseOnly = TRUE){
  # Check that inputs are integer values
  if(!all(sapply(c(ki, kj, M), function(x) x == as.integer(x)))){
    stop("Not all inputs are integer values.")}
  lij <- 1 - exp( -(ki * kj) / (2 * M) )
  return(lij)
}

#' Compute the number of expected connections between a set of vertices.
#'
#' @param v vector of verticies present in graph G.
#' @param G igraph object.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

# Change comp_lambda to depend on verticies, then change method to compute probability based on all combinations of degrees (combn(..., simplify = FALSE)). These changes will need to permiate through to other computations of lambda (bipartite, ...) through bipartite will be only combinations between groups, not within groups. This method is supported by the STRINGdb as well as the Itzkovitz 2003 paper, but may be interpreted differently from the Pradines paper. Need to read more carefully.
# 
comp_lambda <- function(v, G){
  # Check inputs
  stopifnot(class(G) == "igraph")

  # Check for vertices in graph G
  vo <- v[!v %in% V(G)$name]
  if(length(vo) > 0){
    message("Vertices not present in graph: ", paste(vo, collapse = ", "), ".")
    v <- v[v %in% G]
  }

  # Compute lambda
  deg <- igraph::degree(G, v)
  M <- ecount(G)
  if(length(deg) > 0){
    lam <- sum(sapply(combn(deg, 2, simplify = FALSE), function(k, M){
      comp_pij(k[1], k[2], M)}, M = M))
  }else{
    lam <- NA
  }
  return(lam)
}

#' Compute the number of expected connections between two sets of vertices.
#'
#' @param v1 vector of vertex identifiers for list 1.
#' @param v2 vector of vertex identifiers for list 2.
#' @param G igraph of which vertices from v1 and v2 are present.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

comp_lambda_bipartite <- function(v1, v2, G){
  # Check inputs
  stopifnot(class(G) == "igraph")

  # Intersecting vertices in lists need to be removed
  vn <- v1[v1 %in% v2]
  if(length(vn) > 0){
    message(
      "Vertices identified between lists: ", paste(vn, collapse = ", "), ".")
    message("Intersecting vertices removed from lists.")
    v1 <- v1[!v1 %in% vn]
    v2 <- v2[!v2 %in% vn]
  }

  # Check for vertices in graph G
  vu <- unique(c(v1, v2))
  vo <- vu[!vu %in% V(G)$name]
  if(length(vo) > 0){
    message("Vertices not present in graph: ", paste(vo, collapse = ", "), ".")
    v1 <- v1[!v1 %in% vo]
    v2 <- v2[!v2 %in% vo]
  }

  if(any(sapply(list(v1, v2), length) == 0)){
    lam <- NA
  }else{
    # Isolate bipartite edges
    vu <- unique(c(v1, v2))
    sG <- induced_subgraph(G, vu)
    sGel <- as.data.frame(get.edgelist(sG))
    names(sGel) <- c("vi", "vj")
    sGel$vi_ori <- ifelse(sGel$vi %in% v1, "v1", "v2")
    sGel$vj_ori <- ifelse(sGel$vj %in% v1, "v1", "v2")
    sGelb <- sGel[sGel$vi_ori != sGel$vj_ori,]

    # Compute lambda from edge list
    ki <- igraph::degree(G, sGel$vi)
    kj <- igraph::degree(G, sGel$vj)
    M <- ecount(G)
    if(length(ki) > 0 & length(kj) > 0){
      lam <- sum(
        sapply(ki, function(ki, M){
          sapply(kj, function(kj){
            comp_pij(ki, kj, M)})}, M = M))
    }else{
      lam <- NA
    }
    return(lam)
  }
}

#' Compute probability for sparse graph edges.
#'
#' @param lambda expected connections between verticies, computed by 
#' `comp_lambda` or `comp_lambda_bipartite`.
#' @param edgeCnt count of edges present between the verticies tested.
#' @param maxEdgeCnt maximal count of edges possible in graph or subgraph.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

comp_sparse_prob <- function(lambda, edgeCnt, maxEdgeCnt, exact = TRUE){
  if(exact){
    alp <- exp(-lambda) * sum(sapply(0:maxEdgeCnt, function(h, lam){
      exp( h * log(lam) - lfactorial(h))}, lam = lambda))
    prb <- exp(-lambda) * sum(sapply(edgeCnt:maxEdgeCnt, function(h, lam){
      exp( h * log(lam) - lfactorial(h))}, lam = lambda))
  }else{
    prb <- (
      exp(ppois(maxEdgeCnt, lambda, log.p = TRUE)) -
        exp(ppois(edgeCnt, lambda, log.p = TRUE))) /
      (exp(ppois(maxEdgeCnt, lambda, log.p = TRUE)) -
         exp(ppois(0, lambda, log.p = TRUE))
      )
  }
  return(prb)
}

#' Compute edgeset probability.
#'
#' @param v vector of vertex identifiers. Vertices must be present in graph G.
#' @param G igraph object representing the graph and including vertices from v.
#' @param maxEdgeCnt numeric maximum number of edges that could be obsersered
#' given vertices within edgelist.
#' @param exact logical. If TRUE, the approximation of the probability is
#' calculated rather than approximated by `ppois()`.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

comp_edgeset_prob <- function(v, G, exact = TRUE){
  # Check inputs
  stopifnot(class(G) == "igraph")
  v <- unique(v)

  # Calculate lamda, alpha, and approximate the probability of the edge set
  lam <- comp_lambda(v, G)
  z <- ecount(induced_subgraph(G, v))
  M <- length(v)*(length(v)-1)/2
  comp_sparse_prob(lam, z, M, exact)
}

#' Approximate the probability of each vertex within a list belonging to the
#' list.
#'
#' @param v vector of vertex identifiers. Vertices must be present in graph G.
#' @param G igraph object representing the graph and including vertices from v.
#' @param exact logical. If TRUE, the approximation of the probability is
#' calculated rather than approximated by `ppois()`.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

comp_vertex_prob <- function(v, G, exact = TRUE){
  # Check for vertices in graph
  vo <- v[!v %in% V(G)$name]
  if(length(vo) > 0){
    message("Vertices not present in graph: ", paste(vo, collapse = ", "), ".")
    v <- v[!v %in% vo]
  }

  sapply(v, function(vq, vset, G, exact){
    zb <- as.data.frame(get.edgelist(induced_subgraph(G, vset)))
    zb <- nrow(zb[zb[,1] == vq | zb[,2] == vq,])
    vset <- vset[vset != vq]
    lam <- comp_lambda_bipartite(vq, vset, G)
    M <- length(vset) # Since vq has already been removed)
    comp_sparse_prob(lam, zb, M, exact)
  }, vset = v, el = sGel, G = G, exact)
}

#' Approximate the probability of a vertex not within a list belonging to the
#' list.
#'
#' @param v1 vector of vertex identifiers. Vertices must be present in graph G.
#' @param v2 vector of vertex identifiers not in `v1`. These vertices will be
#' tested as a vertex in `v1`.
#' @param G igraph object representing the graph and including vertices from
#' `v1` and `v2`.
#' @param exact logical. If TRUE, the approximation of the probability is
#' calculated rather than approximated by `ppois()`.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

comp_outside_vertex_prob <- function(v1, v2, G, exact = TRUE){
  # All vertices in v1 and v2 must be present in graph G
  vu <- unique(c(v1, v2))
  vo <- vu[!vu %in% names(as.list(V(G)))]
  if(length(vo) > 0){
    message(paste0(
      "Vertices not in supplied graph: ",
      paste(v_not_present, collapse = ", "), "."))
    v1 <- v1[!v1 %in% v_not_present]
    v2 <- v2[!v2 %in% v_not_present]
  }

  # Lists cannot intersect
  # Intersecting object in lists will be removed and messaged.
  vn <- v1[v1 %in% v2]
  if(length(vn) > 0){
    message(
      "Vertices identified between lists: ", paste(vn, collapse = ", "), ".")
    v1 <- v1[!v1 %in% vn]
    v2 <- v2[!v2 %in% vn]
  }

  if(length(v1) == 0 | length(v2) == 0){
    return(NA)
  }else{
    return(sapply(v2, function(vq, vset, G, exact){
      zb <- as.data.frame(get.edgelist(induced_subgraph(G, c(vset, vq))))
      zb <- nrow(zb[zb[,1] == vq | zb[,2] == vq,])
      lam <- comp_lambda_bipartite(vq, vset, G)
      M <- length(vset)
      comp_sparse_prob(lam, zb, M, exact)
    }, vset = v1, G = G, exact = exact))
  }
}

#' Approximate the probability of interactions between two subgraphs, a
#' bipartite comparison.
#'
#' @param v1 vector of vertex identifiers, indicating the vertices present in
#' subgraph 1. Vertex identifiers need to be the names of vertices in G. "v1"
#' vertices will be used to create g1 subgraph.
#' @param v2 vector of vertex identifiers. Same criteria as L1. L1 and L2 for
#' comparison should not have zero intersection. Vertices in both L1 and L2 will
#' be removed during analysis. "v2" vertices will be used to create g2 subgraph.
#' @param G igraph object of which g1 and g2 subgraphs will be created from.
#' @param exact logical. If TRUE, the approximation of the probability is
#' calculated rather than approximated by `ppois()`.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

comp_bipartite_prob <- function(v1, v2, G, exact = TRUE){
  # All vertices in v1 and v2 must be present in graph G
  vu <- unique(c(v1, v2))
  v_not_present <- vu[!vu %in% names(as.list(V(G)))]
  if(length(v_not_present) > 0){
    message(paste0(
      "Vertices not in supplied graph: ",
      paste(v_not_present, collapse = ", "), "."))
    message("Vertices will be removed from input lists.")
    v1 <- v1[!v1 %in% v_not_present]
    v2 <- v2[!v2 %in% v_not_present]
  }

  # Lists cannot intersect
  # Intersecting object in lists will be removed and messaged.
  vn <- v1[v1 %in% v2]
  if(length(vn) > 0){
    message(
      "Vertices identified between lists: ", paste(vn, collapse = ", "), ".")
    message("Intersecting vertices removed from lists.")
    v1 <- v1[!v1 %in% vn]
    v2 <- v2[!v2 %in% vn]
  }

  if(length(v1) == 0 | length(v2) == 0){
    return(NA)
  }else{
    lam <- comp_lambda_bipartite(v1, v2, G)

    # Determine bipartite edges by analyzing subgraph and edgelist
    vu <- unique(c(v1, v2))
    sG <- induced_subgraph(G, vu)
    sGel <- as.data.frame(get.edgelist(sG))
    names(sGel) <- c("vi", "vj")
    sGel$vi_ori <- ifelse(sGel$vi %in% v1, "v1", "v2")
    sGel$vj_ori <- ifelse(sGel$vj %in% v1, "v1", "v2")
    sGelb <- sGel[sGel$vi_ori != sGel$vj_ori,]
    zb <- nrow(sGelb)
    M <- length(unique(v1)) * length(unique(v2))

    if(zb > 0){
      return(comp_sparse_prob(lam, zb, M, exact))
    }else{
      return(NA)
    }
  }
}

#' Calculate the enrichment within an ordered list of verticies.
#'
#' @param v vector of vertex identifiers. Vertices must be present in graph G.
#' @param G igraph object representing the graph and including vertices from v.
#' @param d delta or change distance for which to look for enrichment. Default
#' 20.
#' @param f forward threshold for reference enrichment. Default 100.
#' @param e edge window to consider at the beginning. If not specified, the
#' default value will be equal to the delta `d`.
#' @param exact logical. If TRUE, the approximation of the probability is
#' calculated rather than approximated by `ppois()`.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

comp_enrichment <- function(v, G, d = 20, f = 100, e = NULL, exact = TRUE){
  stopifnot(all(v %in% V(G)$name))
  if(is.null(e)) e <- d
  lims <- seq(e, length(v), d)
  lims <- lims[lims >= e]
  denv <- new.env()
  denv$enrich <- numeric()
  denv$enrichExt <- numeric()

  null <- lapply(lims, function(i, v, G, d, f, e, denv, exact){
    Gec <- ecount(G)
    b <- v[1:max(e, i-d)]
    s <- v[(i-d+1):i]
    if(all(s %in% b)){
      enrich <- comp_edgeset_prob(b, G, exact)
      enrichExt <- enrich
    }else{
      sG <- induced_subgraph(G, s)
      if(ecount(sG) > 0){
        lambdaS <- comp_lambda(as.data.frame(get.edgelist(sG)), G)
        lambdaS <- ifelse(is.na(lambdaS), 0, lambdaS)
        lambdaB <- comp_lambda_bipartite(s, b, G)
        lambdaB <- ifelse(is.na(lambdaB), 0, lambdaB)
        lambdaE <- lambdaB + lambdaS
        edgeCnt <- count_edges(G, s) + count_edges(G, s, b)
        maxEdgeCnt <- length(s)*(length(s)-1)/2 + length(s)*length(b)
        enrich <- comp_sparse_prob(lambdaS, edgeCnt, maxEdgeCnt, exact)

        if(f > 0){
          a <- v[(i+1):min((i+f), length(v))]
          lambdaA <- comp_lambda_bipartite(s, a, G)
          lambdaA <- ifelse(is.na(lambdaA), 0, lambdaA)
          lambdaE <- lambdaE + lambdaA
          edgeCnt <- edgeCnt + count_edges(G, s, a)
          maxEdgeCnt <- maxEdgeCnt + length(s)*length(a)
        }
        enrichExt <- comp_sparse_prob(lambdaE, edgeCnt, maxEdgeCnt, exact)
      }else{
        enrich <- NA
        enrichExt <- NA
      }
    }
    denv$enrich <- c(denv$enrich, enrich)
    denv$enrichExt <- c(denv$enrichExt, enrichExt)
  },
  v = v, G = G, d = d, f = f, e = e, denv = denv, exact = exact)

  return(list(
    "limits" = lims,
    "enrichment" = denv$enrich,
    "enrichmentExt" = denv$enrichExt))
}
