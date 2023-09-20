#' The function markerGSEA calculates gene (marker) set enrichment of
#' cluster-specific variable markers within  bone cell type- or bone
#' tissue-specific marker database.
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
#' @return data.frame with following variables: pathway (name of list variable
#'  with markers used for GSEA) / pval (p-value) / padj (adjusted p-value) /
#'  log2err (log2(error)) / ES (enrichment score) / NES (normalized enrichment
#'  score - useful for multiple comparisons) / size (number of overlaps) /
#'  leadingEdge (overlapping markers) / cluster (cluster number)
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

#' The function plotMarkerGSEA plots results from \code{markerGSEA} function
#
#' @param enrichRes An output from \code{markerGSEA} function in form of a
#'  data.frame.
#' @param topN An optional parameter specifying the number of the most
#'  significant results (based on p-value and NES) to be plotted per cluster.
#' @param plotNES and optional logical variable if bar height should be NES
#'  instead of a -log10(padj), since padj can be really low
#' @return A ggplot object
#'
#' @export
#' @examples
#' varMarkerFile = system.file("extdata", "cluster_variable_features.csv",
#'                             package = "BoneCellType")
#' varMarkers = read.csv(varMarkerFile)
#' enrichRes = markerGSEA(varMarkers)
#' \dontrun{
#' plotMarkerGSEA(enrichRes)}
plotMarkerGSEA = function(enrichRes,
                          topN=Inf,
                          plotNES=FALSE){

  # create lable for plotting (clster and pathway label) -> get negative log10
  # of padj -> arrange data by clusters and significance
  enrichRes = enrichRes %>%
    mutate(label = paste(cluster,pathway, sep = "__")) %>%
    mutate(neg_log10_padj = -log10(padj)) %>%
    group_by(cluster) %>%
    arrange(padj, desc(NES),.by_group = TRUE)

  # factorize labels to be sorted based on cluster number and adjusted
  # p-value of enrichment
  enrichRes$label = factor(enrichRes$label, levels = enrichRes$label)

  # if topN was selected - filter out the top topN results based on padj and
  # NES score
  if (is.finite(topN)) {
    enrichRes = enrichRes %>%
      group_by(cluster) %>%
      arrange(padj, desc(NES),.by_group = TRUE) %>%
      slice(1:topN)
  }


  # should p-values be plotted instead
  if (plotNES) {
    p = ggplot(enrichRes, aes(x = label,
                              y = NES,
                              fill = neg_log10_p)) +
      ylab("NES")
  } else {
    p = ggplot(enrichRes, aes(x = label,
                              y = neg_log10_padj,
                              fill = NES)) +
      ylab(bquote(-log[10](padj)))
  }
  p = p +
    geom_bar(stat= "identity") +
    geom_hline(yintercept = -log10(0.05), color = "grey80", linetype = "dashed") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    xlab("")
  return(p)
}

