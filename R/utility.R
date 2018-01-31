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

#' Assign consistant ID to variable given various alias IDs from a reference.
#'
#' @usage
#'   alias_arbiter(IDs, RefIDs, aliasIDs)
#'   
#'   alias_arbiter(
#'     IDs, RefIDs, aliasIDs, outputIDs, sep = NULL,
#'     remove_absent_IDs = FALSE, quiet = FALSE)
#'
#' @description Assign a reference ID to a vector of IDs that may have alias
#' values. An example, if ID values of "A", "B", "C" were given to three
#' bacterial species, but "B" was later reclassified as "C", then "B" becomes
#' a depreciated alias of "C". This function takes this information and
#' changes the original vector of "A", "B", and "C", to "A", "C", and "C".
#' Aliases can be supplied as a list of character or numeric vectors or as a
#' character vector to split on the `sep` value. IDs need to match RefIDs, but
#' outputIDs can be another ID type, as long as they are indexed in the same
#' order as RefIDs. Any ambiguous aliases will be removed, but will be messaged
#' to the output. Any IDs not matching RefIDs or aliasIDs will be returned as
#' "NA", unless "NULL" is provided to the option `remove_absent_IDs`. For this
#' option, the default value is `FALSE` which will return `NA` for unmatched IDs
#' and `TRUE` will remove all `NA` values. Caution, this produces an output that
#' is not the same length as the input if there are any unmatched or ambiguous
#' IDs. The other option is `NULL` to be supplied, which fills the `NA` values
#' with their original input, even if the alias was ambiguous.
#'
#' This function creates a directed graph where all reference IDs are source
#' nodes and all aliases are children from that source node. Alias IDs are
#' screened for reference IDs and ambiguous IDs first, therefore ensuring the
#' graph is a collection of small subgraphs with single sources and multiple
#' sinks.
#'
#' @param IDs character vector of IDs which need to be consistantly assigned.
#' Must be same type as RefIDs.
#' @param RefIDs character vector of IDs for reference, must be the same length
#' and order as aliasIDs and, if included, outputIDs.
#' @param aliasIDs list of character vectors giving aliases for given RefIDs, or
#' character vector to be split for aliases (must specify the `sep` argument
#' with delimiter).
#' @param outputIDs character or numeric vector for output. Indexed in the same
#' order as the RefIDs and aliasIDs. If not given, RefIDs will be used as
#' the output IDs
#' @param sep character delimiter for aliasIDs if supplied in a character
#' vector. Default is `NULL`.
#' @param remove_absent_IDs logical. Default behavior will change IDs not found
#' in RefIDs or aliases to `NA`. If set to TRUE, IDs not matching RefIDs
#' or aliasIDs will be removed. Alternatively, if set to NULL, IDs not matching
#' RefIDs or aliasIDs will remain as input ID values.
#' @param quiet logical. If TRUE, message output will be silenced. Messages
#' include notification of ambiguous aliases or aliases which belong to two or
#' more reference IDs.
#'
#' @example
#' # Using letters as IDs
#' IDs <- LETTERS[1:10]
#' RefIDs <- c("A", "C", "E", "J")
#' aliasIDs <- c("A,B", "C,D", "E,F,G,H,I", "J")
#' alias_arbiter(IDs, RefIDs, aliasIDs, sep = ",")
#'
#' @author Christopher Nobles, Ph.D.
#' @export

alias_arbiter <- function(IDs, RefIDs, aliasIDs, outputIDs = NULL,
                          sep = NULL, remove_absent_IDs = FALSE, quiet = FALSE){
  # Require packages
  stopifnot(require("igraph"))

  # Check inputs and assign defaults
  stopifnot(any(IDs %in% RefIDs))
  stopifnot(length(RefIDs) == length(aliasIDs))
  if(!is.null(outputIDs)){
    stopifnot(length(outputIDs) == length(RefIDs))
  }else{
    outputIDs <- RefIDs
  }
  if(class(aliasIDs) != "list"){
    stopifnot(any(sapply(aliasIDs, function(x) grepl(sep, x))))
    aliasIDs <- strsplit(aliasIDs, sep, fixed = TRUE)
  }

  # Curate aliases to remove source
  aliasID_df <- data.frame(
    ref = rep(RefIDs, sapply(aliasIDs, length)),
    alias = unlist(aliasIDs),
    stringsAsFactors = FALSE)
  aliasID_df <- aliasID_df[aliasID_df$ref != aliasID_df$alias,]
  aliasID_df <- aggregate(alias ~ ref, data = aliasID_df, function(x){
    paste(x[!x %in% RefIDs], collapse = sep)})
  aliasIDs <- strsplit(aliasID_df$alias, sep)
  names(aliasIDs) <- aliasID_df$ref
  
  # Check for ambiguous aliases
  ambi_alias <- table(unlist(aliasIDs))
  ambi_alias <- names(ambi_alias[ambi_alias > 1])
  if(length(ambi_alias) > 0 & !quiet){
    message(
      "Ambiguous aliases removed : ", paste(ambi_alias, collapse = ", "), ".")
  }

  # Remove ambiguous aliases from aliasIDs
  aliasID_df2 <- data.frame(
    ref = rep(names(aliasIDs), sapply(aliasIDs, length)),
    alias = unlist(aliasIDs, use.names = FALSE),
    stringsAsFactors = FALSE)
  aliasID_df2 <- aliasID_df2[!aliasID_df2$alias %in% ambi_alias,]
  aliasIDs <- split(aliasID_df2$alias, aliasID_df2$ref)

  # Create graph of aliases and reference IDs
  E <- data.frame(
    vi = rep(names(aliasIDs), sapply(aliasIDs, length)),
    vj = unlist(aliasIDs, use.names = FALSE),
    stringsAsFactors = FALSE)
  V <- data.frame(v = unique(c(RefIDs, unlist(aliasIDs))))
  G <- construct_graph(E, V, mode = "directed")

  # Assign IDs to clusters (maybe this should be sources instead of clusters)
  mem <- igraph::clusters(G)$membership
  IDmem <- mem[as.character(IDs)]
  Refmem <- mem[as.character(RefIDs)]

  # Assign output IDs to input IDs and return
  outIDs <- outputIDs[match(IDmem, Refmem)]
  if(is.null(remove_absent_IDs)){
    outIDs <- ifelse(is.na(outIDs), IDs, outIDs)
  }else if(remove_absent_IDs == TRUE){
    outIDs <- outIDs[!is.na(outIDs)]
  }
  return(outIDs)
}
