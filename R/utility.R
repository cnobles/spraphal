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

