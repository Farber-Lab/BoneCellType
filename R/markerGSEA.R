#' The function markerEnrich calculates enrichment of cluster-specific
#' variable markers within  bone cell type- or bone tissue-specific marker
#' database.
#
#' @param varMarkers A data.frame containing variable markers for each cluster
#'  as reported by \code{FindAllMarkers} function from Seurat package.
#' @param markers Optional user provided named list with markers.
#' @param systemLevel Optional string specifying on which systematic level to
#'  perform enrichment, options: "cellType" or "tissue". Default: "cellType".
#'  Use only when not providing your own markers variable.
#' @param byDatabase Optional logical variable specifying if enrichment should
#'  be performed separately for each database (i.e. source  from which markers
#'  were obtained): option TRUE (default). Or if the database column in
#'  \code{markerDatabase} should be ignored and enrichment is done on all
#'  available markers for a given cellType or tissue: option FALSE. Use only
#'  when not providing your own markers variable.
#' @param minSize Optional number specifying a minimal number of markers within
#'  a marker set to be considered (defaults to 5).
#' @param maxSize Optional umber specifying a maximal number of markers within a
#'  marker set to be considered (defaults to 50000).
#' @return data.frame
#'
#' @export
#' @examples
#' varMarkerFile = system.file("extdata", "cluster_variable_features.csv",
#'                             package = "BoneCellType")
#' varMarkers = read.csv(varMarkerFile)
#' GSEAres = markerGSEA(varMarkers)
#'
#' # use own markers
#' markerList = cellType_list
#' GSEAresCustom = markerGSEA(varMarkers, markers=cellType_list)
markerGSEA = function(varMarkers,
                      markers=NULL,
                      systemLevel="cellType",
                      byDatabase=TRUE,
                      minSize = 5,
                      maxSize = 50000){
  # load marked database based on input criteria
  if (is.null(markers)) {
    if (systemLevel == "cellType"){

      if (byDatabase){
        markers = cellTypeSource_list
      } else {
        markers = cellType_list
      }

    } else if (systemLevel == "tissue"){

      if (byDatabase){
        markers = tissueSource_list
      } else {
        markers = tissue_list
      }

    } else {
      stop("Wrong input: systemLevel must be either cellType or tissue.")
    }
  } else if (!is.list(markers)){
    stop("Wrong input: markers must be a list.")
  }


  # go through clusters in varMarkers and for each run GSEA
  for (i in unique(varMarkers$cluster)){

    # filter out given cluster and arrange by log2FC
    # if log2FC is Inf or -Inf set it to 1000 or -1000 respectively (just an
    # arbitrarily high number)
    # also make sure that log2FC is numeric
    testClust = varMarkers %>%
      filter(cluster == i) %>%
      mutate(avg_log2FC = ifelse(is.infinite(avg_log2FC) & avg_log2FC > 0, 1000,
                                 ifelse(is.infinite(avg_log2FC) & avg_log2FC < 0,
                                        -1000, avg_log2FC))) %>%
      arrange(-avg_log2FC)

    # create named vector with log2FC values
    rankedVector = testClust$avg_log2FC
    names(rankedVector) = testClust$gene

    # run fgsea
    fgseaRes = suppressWarnings(fgsea(pathways = markers,
                                      stats    = rankedVector,
                                      minSize  = minSize,
                                      maxSize  = maxSize)) %>%
      mutate(cluster = i)


    if (i == unique(varMarkers$cluster)[1]){
      results = fgseaRes
    } else {
      results = rbind(results, fgseaRes)
    }

  }

  # unlist the genes into a comma separated string
  leadingEdge = c()
  for (i in seq(results$leadingEdge)){
    LE = paste(unlist(results$leadingEdge[[i]]), collapse = ",")
    leadingEdge = c(leadingEdge, LE)
  }
  results$leadingEdge = leadingEdge

  return(results)

}
