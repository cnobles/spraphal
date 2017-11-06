#' Plot subgraph using ggnet2
#'
#' @param v vector of vertex identifiers. All vertex identifiers need to be
#' present in the graph G.
#' @param G igraph object of which v vertices will be used to create a subgraph.
#' @param ... arguments to be passed into ggnet2.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

plot_subgraph <- function(v, G, ...){
  sG.netw <- intergraph::asNetwork(induced_subgraph(G, vids = v))
  ggnet2(net = sG.netw, ...)
}

#' Plot cluster with probablistic coloring of the nodes with relation to the
#' cluster.
#'
#' @param v vector of vertex identifiers. All vertex identifiers need to be
#' present in the graph G.
#' @param G igraph object of which v vertices will be used to create a subgraph.
#' @param ... arguments to be passed into ggnet2.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

plot_cluster <- function(v, G, ...){
  probs <- enrich_vertex(v, G)
  prob.lvls <- c(1, 0.05, 0.01, 0.001, 0)
  probs.cut <- cut(probs, prob.lvls)
  levels(probs.cut) <- c("< 0.001", "< 0.01", "< 0.05", "< 1.0")
  sG.netw <- intergraph::asNetwork(induced_subgraph(G, vids = v))
  set.vertex.attribute(sG.netw, "vertex.prob", as.character(probs.cut))
  #pal <- c("VVS" = "purple", "VS" = "blue", "S" = "green", "NS" = "grey50")
  ggnet2(
    net = sG.netw,
    color = "vertex.prob",
    palette = "Set1",
    legend.position = "right",
    ...) +
    labs(
      title = paste0(
        "Vertices : ", length(v),
        "\nEdges : ", ecount(induced_subgraph(G, vids = v)),
        "\np-value : ",
        format(enrich_edgeset(v, G), digits = 3, scientific = TRUE))) +
    theme(
      panel.background = element_rect(color = "grey"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10))
}

#' Plot bipartite graph using ggnet2
#'
#' @param v1 vector of vertex identifiers for subgraph g1. All vertex
#' identifiers need to be present in the graph G.
#' @param v2 vector of vertex identifiers for subgraph g2. All vertex
#' identifiers need to be present in the graph G.
#' @param G igraph object of which v1 and v2 vertices will be used to create a
#' subgraph.
#' @param ... arguments to be passed into ggnet2.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

plot_bipartite <- function(v1, v2, G, ..., quiet = FALSE){
  vn <- v1[v1 %in% v2]
  if(length(vn) > 0 & !quiet){
    message("Removing intersecting vertices: ", paste(vn, collapse = ", "), ".")
  }

  v1 <- v1[!v1 %in% vn]
  v2 <- v2[!v2 %in% vn]

  sGel <- as.data.frame(get.edgelist(induced_subgraph(G, c(v1, v2))))
  names(sGel) <- c("vi", "vj")
  sGel$vi_ori <- ifelse(sGel$vi %in% v1, "v1", "v2")
  sGel$vj_ori <- ifelse(sGel$vj %in% v1, "v1", "v2")
  sGel <- sGel[sGel$vi_ori != sGel$vj_ori,]

  df <- data.frame(
    vu = c(v1, v2), grp = c(rep("A", length(v1)), rep("B", length(v2))))

  sG.netw <- intergraph::asNetwork(induced_subgraph(G, vids = c(v1, v2)))
  set.vertex.attribute(
    sG.netw, "Group", df$grp[
      match(get.vertex.attribute(sG.netw, "vertex.names"), df$vu)])

  ggnet2(net = sG.netw, color = "Group", palette = "Set2", ...) +
    labs(title = paste0("Vertices: A = ", length(v1), " B = ", length(v2), "\n",
                        "Edges: A = ", ecount(induced_subgraph(G, v1)),
                        " B = ", ecount(induced_subgraph(G, v2)), "\n",
                        "Bipartite edges: ", nrow(sGel), "\n",
                        "Interaction p-value: ", format(
                          enrich_bipartite(v1, v2, G), digits = 3,
                          scientific = TRUE))) +
    theme(
      panel.background = element_rect(color = "grey"),
      title = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10))
}

#' Plot graph with a list of vertex lists using ggnet2
#'
#' @param Lc list of character vectors holding information on vertex ids for
#' different clusters or groups of vertices. Each object in the list will be
#' given a different color. Additionally, names of the list will be used in the
#' legend, unless no names are supplied. If a vertex in named on more than one
#' list, it will be colored black and pushed into an intersecting group.
#' @param G igraph object from which vertices in Lc will be used to create a
#' subgraph.
#' @param set.v character, Set name to be passed on to RColorBrewer::brewer.pal
#' to select color palette for vertices. If name supplied is not in
#' brewer.pal.info, then the name will be treated as a color and passed as a
#' constant color for vertices.
#' @param set.e character, Set name to be passed on to RColorBrewer::brewer.pal
#' to select color palette for edges. If name supplied is not in
#' brewer.pal.info, then the name will be treated as a color and passed as a
#' constant color for edges.
#' @param ... arguments to be passed into ggnet2.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

plot_clusters <- function(Lc, G, set.v = "Set1", set.e = "grey50", ...){
  # Required R-packages
  packs <- c("igraph", "sna", "intergraph", "ggnet", "spraphal", "RColorBrewer",
             "ggplot2", "scales")
  stopifnot(any(sapply(packs, require, character.only = TRUE)))

  # Lc should be a list of vertex vectors, create list of 1 if given a vector.
  if(class(Lc) != "list"){
    stopifnot(is.vector(Lc))
    Lc <- list(unique(Lc))
  }else{
    Lc <- lapply(Lc, unique)
  }

  # Lc should have names, if not a component of the list, then give letters.
  if(is.null(names(Lc))) names(Lc) <- LETTERS[1:length(Lc)]
  L_names <- names(Lc)

  # Input list needs to contain vertices within G
  vu <- unique(unlist(Lc))
  vm <- vu[!vu %in% names(V(G))]
  if(length(vm) > 0){
    message(
      "Removing missing vertices from input list: ",
      paste(vm, collapse = ", "), ".")
    Lc <- lapply(Lc, function(x) x[!x %in% vm])
    names(Lc) <- L_names
  }
  if(sum(sapply(Lc, length)) == 0){
    stop("No remaining vertices within the vertex list. Vertex names likely
         do not match vertex names in graph G.")
  }

  # Label intersecting vertices in own group.
  vn <- table(unlist(Lc))
  vn <- names(vn)[vn > 1]
  Lc <- lapply(Lc, function(x) x[!x %in% vn])
  Lc <- c(Lc, list("Vn" = vn))

  # Generate subgraph and edge list
  sG <- induced_subgraph(G, unlist(Lc))
  sGel <- as.data.frame(get.edgelist(sG))
  names(sGel) <- c("vi", "vj")
  sGel$vi_ori <- sapply(sGel$vi, function(x) names(Lc)[grep(x, Lc)])
  sGel$vj_ori <- sapply(sGel$vj, function(x) names(Lc)[grep(x, Lc)])

  # Annotate colors to edges and vertices
  sGel$edge.type <- ifelse(
    with(sGel, vi_ori != vj_ori & vi_ori != "Vn" & vj_ori != "Vn"),
    "bipartite",
    ifelse(
      with(sGel, vi_ori == vj_ori & vi_ori != "Vn" & vj_ori != "Vn"),
      "interedge",
      "shared"))
  if(set.e %in% row.names(brewer.pal.info)){
    set.e.max <- brewer.pal.info$maxcolors[
      match(set.e, row.names(brewer.pal.info))]
    edge_colors <- colorRampPalette(brewer.pal(set.e.max, set.e))(
      length(unique(sGel$edge.type)))
  }else{
    edge_colors <- rep(set.e, length(unique(sGel$edge.type)))
  }
  names(edge_colors) <- unique(sGel$edge.type)

  if(set.v %in% row.names(brewer.pal.info)){
    set.v.max <- brewer.pal.info$maxcolors[
      match(set.v, row.names(brewer.pal.info))]
    Lc_colors <- colorRampPalette(brewer.pal(set.v.max, set.v))(
      length(Lc)-1)
  }else{
    Lc_colors <- rep(set.v, length(Lc)-1)
  }
  Lc_colors <- c(Lc_colors, "#000000")
  names(Lc_colors) <- names(Lc)

  sGvl <- data.frame(
    "vi" = unlist(Lc),
    "grp" = rep(names(Lc), sapply(Lc, length)),
    row.names = NULL)
  sGvl <- sGvl[match(V(sG)$name, sGvl$vi),]

  # Apply attributes of vertices and edges to sG
  E(sG)$edge.type <- sGel$edge.type
  E(sG)$edge.color <- edge_colors[sGel$edge.type]
  V(sG)$Group <- as.character(sGvl$grp)
  edge_type_count <- table(sGel$edge.type)
  if(length(edge_type_count) == 3){
    if(!"bipartite" %in% names(edge_type_count)){
      edge_type_count <- c(edge_type_count, "bipartite" = 0)}
    if(!"interedge" %in% names(edge_type_count)){
      edge_type_count <- c(edge_type_count, "interedge" = 0)}
    if(!"shared" %in% names(edge_type_count)){
      edge_type_count <- c(edge_type_count, "shared" = 0)}
  }

  # Construct plot
  if(length(Lc) > 1){
    p <- ggnet2(
        net = sG,
        color = "Group", palette = Lc_colors,
        edge.color = "edge.color", ...)
  }else{
    p <- ggnet2(
      net = sG,
      color = Lc_colors[1],
      edge.color = "edge.color", ...)
  }
  p + labs(
      title = paste0(
        "Vertices: Total = ", nrow(sGvl), " : ", paste(
          paste0(names(Lc), " = ", sapply(Lc, length)), collapse = ", "), "\n",
        "Edges: Total = ", nrow(sGel), " : Bipartite = ",
          edge_type_count["bipartite"], ", Interedge = ",
          edge_type_count["interedge"], ", Shared = ",
          edge_type_count["shared"])) +
    theme(
      panel.background = element_rect(color = "grey"),
      title = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10))
}
