#' Test graph for sparseness
#' 
#' @param G igraph object representing the graph to test.
#' @param threshold numeric, between 0 and 1 indicating the proportion of which
#' to consider much greater than. For example, if 10 is much greater than 1, 
#' then the threshold is considered to be 0.1. Default threshold is 0.1.
#' 
#' @author Christopher Nobles, Ph.D.
#' @export

sparseness <- function(G, threshold = 0.1){
  stopifnot(require("igraph"))
  V <- V(G)$name
  el <- as.data.frame(get.edgelist(G))
  names(el) <- c("vi", "vj")
  el$ki <- igraph::degree(G, v = el$vi)
  el$kj <- igraph::degree(G, v = el$vj)
  el$kikj <- el$ki * el$kj
  
  # Quantify sparseness
  M <- ecount(G)
  maxM <- ifelse(
    igraph::is_directed(G),
    length(V) * (length(V) - 1 ),
    length(V) * (length(V) - 1 ) / 2)
  denG <- igraph::graph.density(G)
  kimax <- max(c(el$ki, el$kj))
  kikjmax <- max(el$kikj)
  pct_above_ki <- 100 * ( 
    sum(as.integer(igraph::degree(G, V) >= threshold * length(V))) ) / length(V)
  pct_above_kikj <- 100 * ( 
    sum(as.integer(el$kikj >= threshold * M)) ) / M
  
  # Ouput Results
  return(list(
    "N" = length(V),
    "M" = M,
    "maxM" = maxM,
    "Density" = denG,
    "kiMax" = kimax,
    "pct_above_ki_threshold" = pct_above_ki,
    "kikjMax" = kikjmax,
    "pct_above_kikj_threshold" = pct_above_kikj))
}

#' Match gene names within data to vertex identifiers with a reference data set.
#' 
#' @param data data.frame with column specified by the `names` parameter.
#' @param names character specifying the column name in data containing the 
#' `names` of to be mapped from.
#' @param refdata data.frame with reference data. Needs to contain the columns
#' referred to by `refcol` and `refout`. 
#' @param refcol character the name of the column with matching information for
#' the `names` parameter.
#' @param refout character the name of the column for in the reference with the 
#' information to be mapped to data.
#' @param outname character the name given to the new information mapped to the 
#' `data` data.frame.
#' @param removeUnmapped logical indicating if rows should be removed if
#' unmapped.
#' @param quiet logical indicating if messages should be printed or not.
#' 
#' @author Christopher Nobles, Ph.D.
#' @export

map_names <- function(data, names = "gene_name", 
                      refdata, refcol = "preferred_name", 
                      refout = "protein_external_id", 
                      outname = "stringID", removeUnmapped = TRUE, 
                      quiet = FALSE){
  # Check input data sets for data.frame structure
  stopifnot(all(sapply(list(data, refdata), class) == "data.frame"))
  
  # Identify which columns to use
  dataColNum <- match(names, names(data))
  refColNum <- match(refcol, names(refdata))
  refOutNum <- match(refout, names(refdata))
  
  # Match names from data to reference
  matched_order <- match(data[,dataColNum], refdata[,refColNum])
  
  # Annotate data with output identifier
  data[,outname] <- refdata[,refout][matched_order]
  
  # Message unmapped identifiers
  if(removeUnmapped & !quiet){
    unmapped <- sum(as.integer(is.na(matched_order)))
    unmappedPct <- round(100 * unmapped / nrow(data), digits = 2)
    message(paste0(
      "Unmapped: ", unmappedPct, "% (", unmapped, " out of ", nrow(data), ")"))
  }
  
  # Remove unmapped rows if requested
  if(removeUnmapped){
    data <- data[!is.na(data[,outname]),]
  }
  
  return(data)
}

#' Construct graph from data.frame
#' 
#' @param df data.frame where the initial two columns represent an edgelist.
#' @param V vector of vertex names / identifiers. Vertices not present in the 
#' edgelist will be included as vertices with no edges in the output graph. If 
#' nothing supplied, only vertices present in edgelist data.frame will be
#' present in the graph.
#' @param dfcols numeric vector of length 2, specifying the columns of the 
#' data.frame that contain the edgelist information.
#' @param mode character. Indicates if the graph is a directed or undirected 
#' 
#' @author Christopher Nobles, Ph.D.
#' @export

construct_graph <- function(df, V = NULL, dfcols = c(1,2), mode = "undirected"){
  required <- c("igraph")
  stopifnot(all(sapply(
    required, require, character.only = TRUE, quietly = TRUE)))
  
  # Determine vertices of graph
  if(is.null(V)){
    V <- unique(c(df[,dfcols[1]], df[,dfcols[2]]))
  }else{
    V <- unique(c(df[,dfcols[1]], df[,dfcols[2]], V))
  }
  
  # Construct initial graph
  df_el <- data.frame(
    vi = factor(df[,dfcols[1]], levels = V),
    vj = factor(df[,dfcols[2]], levels = V))
  
  return(graph_from_adjacency_matrix(adjmatrix = table(df_el), mode = mode))
}

