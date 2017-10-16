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
        "p-value = ", 
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
                        "Interaction p-value: ", enrich_bipartite(v1, v2, G))) +
    theme(
      panel.background = element_rect(color = "grey"),
      title = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10))
}
