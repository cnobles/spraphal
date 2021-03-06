#' Construct graph from data.frame
#'
#' @param E data.frame with edgelist information, by default the first two
#' columns should contain the vertex connection data. This can be changed with
#' the `ecols` parameter. Additional data will be assimilated into the graph as
#' edge attributes.
#' @param V data.frame with vertex information. Vertices not present in the
#' edgelist will be included as vertices with no edges in the output graph. If
#' nothing supplied, only vertices present in edgelist data.frame will be
#' present in the graph. Additional data will be assimilated into the graph as
#' vertex attributes.
#' @param ecols numeric vector of length 2, specifying the columns of the
#' data.frame E that contain the edgelist information.
#' @param vcol numeric vector of length 1, specifying the column with vertex
#' identifier information within the data.frame V.
#' @param mode character. Indicates if the graph is a directed or undirected
#'
#' @author Christopher Nobles, Ph.D.
#' @export

construct_graph <- function(E, V = NULL, ecols = c(1,2), vcol = 1,
                            mode = "undirected"){
  # Check inputs
  stopifnot(require("igraph"))
  stopifnot(class(E) == "data.frame")
  stopifnot(
    all(c(is.numeric(ecols), length(ecols) == 2, max(ecols) <= ncol(E))))
  if(!is.null(V)){
    stopifnot(class(V) == "data.frame")
    stopifnot(all(c(is.numeric(vcol), length(vcol) == 1, max(vcol) <= ncol(V))))
  }
  stopifnot(mode %in% c("directed", "undirected"))

  # Format data.frame E
  if(any(!ecols %in% c(1,2))){
    e_colnames <- names(E)
    E <- cbind(E[,ecols], E[,-ecols])
    names(E) <- c(e_colnames[ecols], e_colnames[-ecols])
  }

  # Format data.frame V
  if(!is.null(V) & vcol != 1){
    v_colnames <- names(V)
    V <- cbind(V[,vcol], V[,-vcol])
    names(V) <- c(v_colnames[vcol], v_colnames[-vcol])
  }

  # Construct graph
  return(graph_from_data_frame(
    d = E, directed = (mode == "directed"), vertices = V))
}

#' Count edges within a graph or edges between two non-overlapping graphs
#'
#' @param G igraph which contains vertex identified by `v1` and `v2`.
#' @param v1 vector of vertex identifiers.
#' @param v2 vector of vertex identifiers. `NULL` by default, but when a second
#' vector is supplied, the fuction will return the number of edges connecting
#' the two subgraphs indicated by the vertex lists.
#' @param edges logical. If TRUE, the edges rather than the count will be
#' returned as a data.frame with two columns named "vi" and "vj".
#'
#' @author Christopher Nobles, Ph.D.
#' @export

count_edges <- function(G, v1, v2 = NULL, edges = FALSE){
  stopifnot(require("igraph"))
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
    message("Noted vertices removed from list.")
    v1 <- v1[!v1 %in% vo]
    v2 <- v2[!v2 %in% vo]
  }

  if(is.null(v2)){
    if(!edges){
      return(ecount(induced_subgraph(G, v1)))
    }else{
      el <- as.data.frame(as_edgelist(induced_subgraph(G, v1)))
      names(el) <- c("vi", "vj")
      return(el)
    }
  }else{
    el <- as.data.frame(as_edgelist(induced_subgraph(G, unique(c(v1, v2)))))
    names(el) <- c("vi", "vj")
    el$vi_ori <- ifelse(el$vi %in% v1, "v1", "v2")
    el$vj_ori <- ifelse(el$vj %in% v1, "v1", "v2")
    el <- el[el$vi_ori != el$vj_ori,]
    if(!edges){
      return(nrow(el))
    }else{
      return(el[,c("vi", "vj")])
    }
  }
}

#' Prune vertices from vertex list based on p-value cutoff / enrichment
#'
#' @param G igraph which contains vertex identifiers present in `v`.
#' @param v vector of vertex identifiers to be considered for pruning.
#' @param cutoff numeric. Cutoff enrichment value (p-value) for which to prune
#' vertices from the list.
#' @param collect logical. If TRUE, the output will contain a listed object with
#' the names of `pruned_list` and `removed_vertices`.
#' @author Christopher Nobles, Ph.D.
#' @export

prune_by_enrich <- function(G, v, cutoff, collect = FALSE){
  stopifnot(require("igraph"))
  stopifnot(class(G) == "igraph")

  # Check for vertices in graph G
  vo <- v[!v %in% V(G)$name]
  if(length(vo) > 0){
    message("Vertices not present in graph: ", paste(vo, collapse = ", "), ".")
    message("Noted vertices removed from list.")
    v <- v[!v %in% vo]
  }

  # Determine the probability of each vertex belonging to the list
  p <- comp_vertex_prob(v, G)
  p <- p[!is.na(p)]
  p <- p[p <= cutoff]
  if(!collect){
    return(names(p))
  }else{
    return(
      list("pruned_list" = names(p), "removed_vertices" = v[!v %in% names(p)]))
  }
}

#' Add vertices to a list based on p-value cutoff / enrichment
#'
#' @param G igraph which contains vertex identifiers present in `v1` and `v2`.
#' @param v1 vector of vertex identifiers. This list will be added to from `v2`
#' if there are any vertices passing the enrichment criteria.
#' @param v2 vector of vertex identifiers to potentially add to `v1`.
#' @param cutoff numeric. Cutoff enrichment value (p-value) for which to add
#' vertices from `v2` to `v1` the list.
#' @param collect logical. If TRUE, the output will contain a listed object with
#' the names of `adjusted_list` and `remaining_vertices`.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

append_by_enrich <- function(G, v1, v2, cutoff, collect = FALSE){
  stopifnot(require("igraph"))
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
    message("Noted vertices removed from list.")
    v1 <- v1[!v1 %in% vo]
    v2 <- v2[!v2 %in% vo]
  }

  # Determine the probability of the vertex associating with the list
  p <- comp_outside_vertex_prob(v1, v2, G)
  p <- p[!is.na(p)]
  p <- p[p <= cutoff]
  if(!collect){
    return(c(v1, names(p)))
  }else{
    return(list(
      "adjusted_list" = c(v1, names(p)),
      "remaining_vertices" = v2[!v2 %in% names(p)]))
  }
}

#' Append vertices in a list to vertex lists depending on enrichment
#'
#' @param G igraph which contains vertex identifiers present in `L` and `v`.
#' @param L list of vectors with vertex identifiers. This list will function as
#' the base lists for which vertices from `v` will be added to based on the best
#' enrichment.
#' @param v vector of vertex identifiers to be added to vertex lists in `L`.
#' Vertices will only be added if they are below the `cutoff` p-value. If there
#' are multiple possibilities, vertices will only be added to the most
#' significant match.
#' @param cutoff numeric. Cutoff enrichment value (p-value) for which to add
#' vertices from `v` to vectors in `L`.
#' @param collect logical. If TRUE, the output will contain a listed object with
#' the names of `adjusted_lists` and `remaining_vertices`.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

append_vertices_to_lists <- function(G, L, v, cutoff, collect = FALSE){
  stopifnot(require("igraph"))
  stopifnot(class(G) == "igraph")

  # Intersecting vertices in lists need to be removed
  vn <- v[v %in% unique(unlist(L))]
  if(length(vn) > 0){
    stop(
      "Vertices within 'v' are already present in 'L'. Intersecting vertices: ",
      paste(vn, collapse = ", "), ".")
  }

  # Check for vertices in graph G
  vu <- unique(c(v, unique(unlist(L))))
  vo <- vu[!vu %in% V(G)$name]
  if(length(vo) > 0){
    stop(
      "Not all vertices present in graph 'G'. Missing vertices: ",
      paste(vo, collapse = ", "), ".")
  }

  # Construct matrix of potential assignment, cols -> L, rows -> v
  # Identify best match and associated pvalue
  assign_matrix <- sapply(L, comp_outside_vertex_prob, v2 = v, G = G)
  # Errors out if the matrix is only NA
  na_rows <- sapply(1:nrow(assign_matrix), function(i){
    sum(as.integer(is.na(assign_matrix[i,])))}) == length(L)
  assign_matrix <- assign_matrix[!na_rows,]
  if(nrow(assign_matrix) > 0){
    best_match <- sapply(1:nrow(assign_matrix), function(i){
      which(assign_matrix[i,] == min(assign_matrix[i,], na.rm = TRUE))})
    min_pvals <- sapply(1:nrow(assign_matrix), function(i){
      min(assign_matrix[i,], na.rm = TRUE)})

    # Filter and select which vertices from v to add to L lists
    v_df <- data.frame(
      v = row.names(assign_matrix), list = best_match, pval = min_pvals)
    v_df <- v_df[v_df$pval <= cutoff,]
    L_adj <- lapply(1:length(L), function(i){
      l <- L[[i]]
      df <- v_df[v_df$list == i,]
      return(c(l, as.character(df$v)))
    })
  }else{
    L_adj <- L
  }

  if(!collect){
    return(L_adj)
  }else{
    return(list(
      "adjusted_lists" = L_adj,
      "remaining_vertices" = v[!v %in% unique(unlist(L_adj))]
    ))
  }
}
