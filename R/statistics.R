#' Compute the probability of vertex i connecting to vertex j
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

#' Compute the number of expected connections between a set of vertices
#'
#' @param el data.frame/matrix edgelist where column 1 and column 2 contain the
#' lists of connected vertices of interest.
#' @param G igraph graph of which vertices from `el` are present.
#' @param elcols numeric vector of length 2, indicating the two columns of the
#' input data.frame or matrix (`el`) which contain the edge list data, defaults
#' to 1 and 2.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

comp_lambda <- function(el, G, elcols = c(1,2)){
  # Check inputs
  stopifnot(class(el) %in% c("data.frame", "matrix"))
  stopifnot(class(G) == "igraph")
  stopifnot(
    all(c(is.numeric(elcols), length(elcols) == 2, max(elcols) <= ncol(el))))

  # Check for vertices in graph G
  vu <- unique(c(el[,elcols[1]], el[,elcols[2]]))
  vo <- vu[!vu %in% V(G)$name]
  if(length(vo) > 0){
    message("Vertices not present in graph: ", paste(vo, collapse = ", "), ".")
    message("Noted vertices removed from edge list.")
    el <- el[el[,!elcols[1]] %in% vo | !el[,elcols[2]] %in% vo,]
  }

  # Compute lambda from edge list
  ki <- igraph::degree(G, el[,elcols[1]])
  kj <- igraph::degree(G, el[,elcols[2]])
  M <- ecount(G)
  l <- sum(mapply(comp_pij, ki = ki, kj = kj, MoreArgs = list(M = M)))
  return(l)
}

#' Compute the number of expected connections between two sets of vertices
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
    message("Noted vertices removed from edge list.")
    v1 <- v1[!v1 %in% vo]
    v2 <- v2[!v2 %in% vo]
  }

  if(any(sapply(list(v1, v2), length) == 0)){
    return(NA)
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
    l <- sum(mapply(comp_pij, ki = ki, kj = kj, MoreArgs = list(M = M)))
    return(l)
  }
}

#' Compute edgeset probability
#'
#' @param el data.frame/matrix edgelist where column 1 and column 2 contain the
#' lists of connected vertices of interest.
#' @param G igraph graph of which vertices from `el` are present.
#' @param maxEdgeCnt numeric maximum number of edges that could be obsersered
#' given vertices within edgelist.
#' @param elcols numeric vector of length 2, indicating the two columns of the
#' input data.frame or matrix (`el`) which contain the edge list data, defaults
#' to 1 and 2.
#' @param exact logical. If TRUE, the approximation of the probability is
#' calculated rather than approximated by `ppois()`.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

comp_edgeset_prob <- function(el, G, maxEdgeCnt,
                              elcols = c(1,2), exact = TRUE){
  # Check inputs
  stopifnot(class(el) %in% c("data.frame", "matrix"))
  stopifnot(class(G) == "igraph")
  stopifnot(
    all(c(is.numeric(elcols), length(elcols) == 2, max(elcols) <= ncol(el))))

  # Check for vertices in graph G
  vu <- unique(c(el[,elcols[1]], el[,elcols[2]]))
  vo <- vu[!vu %in% V(G)$name]
  if(length(vo) > 0){
    message("Vertices not present in graph: ", paste(vo, collapse = ", "), ".")
    message("Noted vertices removed from edge list.")
    el <- el[el[,!elcols[1]] %in% vo | !el[,elcols[2]] %in% vo,]
  }

  # Calculate lamda, alpha, and approximate the probability of the edge set
  vu <- unique(c(el[,elcols[1]], el[,elcols[2]]))
  lam <- comp_lambda(el, G, elcols = elcols)

  if(exact){
    alp <- exp(-lam) * sum(sapply(0:maxEdgeCnt, function(h, lam){
      exp( h * log(lam) - lfactorial(h))}, lam = lam))
    prb <- exp(-lam) * sum(sapply(nrow(el):maxEdgeCnt, function(h, lam){
      exp( h * log(lam) - lfactorial(h))}, lam = lam))
  }else{
    prb <- (
      exp(ppois(maxEdgeCnt, lam, log.p = TRUE)) -
        exp(ppois(nrow(el), lam, log.p = TRUE))) /
      (exp(ppois(maxEdgeCnt, lam, log.p = TRUE)) -
         exp(ppois(0, lam, log.p = TRUE))
      )
  }
  return(prb)
}

#' Approximate the probability of observing an edge set
#'
#' @param v vector of vertex identifiers. Vertices must be present in graph G.
#' @param G igraph object representing the graph and including vertices from v.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

enrich_edgeset <- function(v, G){
  required <- c("igraph", "stats")
  stopifnot(all(sapply(
    required, require, character.only = TRUE, quietly = TRUE)))

  # Check for vertices in graph
  vo <- v[!v %in% names(as.list(V(G)))]
  if(length(vo) > 0){
    message("Vertices not present in graph: ", paste(vo, collapse = ", "), ".")
    v <- v[!v %in% vo]
  }

  # Extract subgraph from G containing all vertices supplied
  sG <- induced_subgraph(G, vids = v)

  # Extract edgelist
  sGel <- as.data.frame(get.edgelist(sG))
  names(sGel) <- c("vi", "vj")

  if(nrow(sGel) == 0){
    return(NA)
  }else{
    prob <- comp_edgeset_prob(
      el = sGel, G = G, maxEdgeCnt = ((length(v) * (length(v) - 1)) / 2))
    return(prob)
  }
}

#' Approximate the probability of each vertex within a list belonging to the
#' list
#'
#' @param v vector of vertex identifiers. Vertices must be present in graph G.
#' @param G igraph object representing the graph and including vertices from v.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

enrich_vertex <- function(v, G){
  required <- c("igraph", "stats")
  stopifnot(all(sapply(
    required, require, character.only = TRUE, quietly = TRUE)))

  # Check for vertices in graph
  vo <- v[!v %in% names(as.list(V(G)))]
  if(length(vo) > 0){
    message("Vertices not present in graph: ", paste(vo, collapse = ", "), ".")
    v <- v[!v %in% vo]
  }

  # Extract subgraph and set up edge list
  sG <- induced_subgraph(G, vids = v)
  sGel <- as.data.frame(get.edgelist(sG))
  names(sGel) <- c("vi", "vj")

  sapply(v, function(vq, vset, el, G){
    # Identify edges only associated with vertex vq
    el <- el[el$vi == vq | el$vj == vq,]
    if(nrow(el) == 0){
      return(NA)
    }else{
      prob <- comp_edgeset_prob(
        el = el, G = G, maxEdgeCnt = (length(vset) - 1))
      return(prob)
    }
  }, vset = v, el = sGel, G = G)
}

#' Approximate the probability of interactions between two subgraphs, a
#' bipartite comparison
#'
#' @param v1 vector of vertex identifiers, indicating the vertices present in
#' subgraph 1. Vertex identifiers need to be the names of vertices in G. "v1"
#' vertices will be used to create g1 subgraph.
#' @param v2 vector of vertex identifiers. Same criteria as L1. L1 and L2 for
#' comparison should not have zero intersection. Vertices in both L1 and L2 will
#' be removed during analysis. "v2" vertices will be used to create g2 subgraph.
#' @param G igraph object of which g1 and g2 subgraphs will be created from.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

enrich_bipartite <- function(v1, v2, G){
  required <- c("igraph", "stats")
  stopifnot(all(sapply(
    required, require, character.only = TRUE, quietly = TRUE)))

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
    return(1)
  }else{
    # Extract subgraph from G containing all vertices from lists
    vu <- unique(c(v1, v2))
    sG <- induced_subgraph(G, vu)

    # Extract edgelist and filter for bipartite edges
    sGel <- as.data.frame(get.edgelist(sG))
    names(sGel) <- c("vi", "vj")
    sGel$vi_ori <- ifelse(sGel$vi %in% v1, "v1", "v2")
    sGel$vj_ori <- ifelse(sGel$vj %in% v1, "v1", "v2")
    sGelb <- sGel[sGel$vi_ori != sGel$vj_ori,]
    v1 <- unique(
      c(sGelb$vi[sGelb$vi_ori == "v1"], sGelb$vj[sGelb$vj_ori == "v1"]))
    v2 <- unique(
      c(sGelb$vi[sGelb$vi_ori == "v2"], sGelb$vj[sGelb$vj_ori == "v2"]))

    if(nrow(sGelb) > 0){
      prob <- comp_edgeset_prob(
        el = sGelb, G = G, maxEdgeCnt = (length(v1) * length(v2)))
      return(prob)
    }else{
      return(NA)
    }
  }
}
