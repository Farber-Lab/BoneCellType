#' The function markerEnrich calculates enrichment of cluster-specific
#' variable markers within  bone cell type- or bone tissue-specific marker
#' database.
#
#' @param varMarkers A data.frame containing variable markers for each cluster
#'  as reported by \code{FindAllMarkers} function from Seurat package.
#' @param markerDatabase A data.frame containing bone-specific markers. The
#'  data.frame should contain following columns: tissue (specifies tissue),
#'  cellType (specifies cell type), gene (contains markers), database (specifies
#'  source from which markers were obtained)
#' @param systemLevel Optional string specifying on which systematic level to
#'  perform enrichment, options: "cellType" or "tissue". Default: "cellType".
#' @param byDatabase Optional logical variable specifying if enrichment should
#'  be performed separately for each database (i.e. source  from which markers
#'  were obtained): option TRUE (default). Or if the database column in
#'  \code{markerDatabase} should be ignored and enrichment is done on all
#'  available markers for a given cellType or tissue: option FALSE.
#' @param topN Optional integer specifying how many top variable markers per
#'  cluster should be used for the enrichment analysis. Default: all variable
#'  markers are used.
#' @param sortBy Optional string used only if topN option is used. Specifies if
#'  log2(fold change): option "log2FC" (default) or adjusted p-value: option
#'  "padj" is used to select the topN variable markers.
#' @return A data.frame with marker enrichment results. Conatins following
#'  columns: cluster (cluster number) / cellType or tissue (based on input
#'  selection) / database (if baDatabase=TRUE) / p (p-value for enrichment) /
#'  oddsRatio (odds ratio for enrichment) / n_olaps (number of overlaps between
#'  variable markers and a given marker set) / overlaps (specific list of
#'  markers found in both variable markers and a given marker set) / padj (
#'  Benjamini Hochberg adjusted p-value) / neg_log10_padj (-log10(padj))
#'
#' @export
#' @examples
#' varMarkerFile = system.file("extdata", "cluster_variable_features.csv",
#'                             package = "BoneCellType")
#' varMarkers = read.csv(varMarkerFile)
#' markerDatabase = markerDatabase
#' enrichRes = markerEnrich(varMarkers, markerDatabase)
markerEnrich = function(varMarkers,
                        markerDatabase,
                        systemLevel="cellType",
                        byDatabase=TRUE,
                        topN=Inf,
                        sortBy="log2FC"){

  # if topN was specified filter topN variable markers based on selected
  # criteria
  if (is.finite(topN)){
    if(sortBy=="padj"){

      varMarkers = varMarkers %>%
        group_by(cluster) %>%
        top_n(n = topN, wt = avg_log2FC)

    } else if (sortBy=="log2FC"){

      varMarkers = varMarkers %>%
        group_by(cluster) %>%
        top_n(n = topN, wt = avg_log2FC)

    } else {
      stop("Wrong input: sortBy must be either padj or log2FC")
    }

  }

  # get marker universe = all unique genes (= markers) in variable markers and
  # bone maker database
  geneUniverse = unique(c(varMarkers$gene, markerDatabase$gene))

  # Create "testTable" from marker table based on systemLevel and byDatabase
  # input. testTable groups markers into subgroups based on these settings
  if (systemLevel == "tissue"){

    # and group by database?
    if (byDatabase){
      testTable = markerDatabase %>%
        select(tissue, database, gene) %>%
        distinct() %>%
        mutate(myGroup = paste(tissue, database, sep = "__"))
    } else {
      testTable = markerDatabase %>%
        select(tissue, gene) %>%
        distinct() %>%
        mutate(myGroup = tissue)
    }

  } else if (systemLevel == "cellType"){

    # and group by database?
    if (byDatabase){
      testTable = markerDatabase %>%
        select(cellType, database, gene) %>%
        distinct() %>%
        mutate(myGroup = paste(cellType, database, sep = "__"))
    } else {
      testTable = markerDatabase  %>%
        select(cellType, gene) %>%
        distinct() %>%
        mutate(myGroup = cellType)
    }

  } else {
    stop("Wrong input: systemLevel must be either cellType or tissue.")
  }

  # calculate enrichment per cluster and per marker group (i.e. per previously
  # created myGroup variable in testTable)
  for (clust in unique(varMarkers$cluster)){

    # extract variable features for a given cluster
    subVar = varMarkers %>%
      filter(cluster == clust)

    for (helpGroup in unique(testTable$myGroup)){

      # extract marker genes for a given group
      subMark = testTable %>%
        filter(myGroup == helpGroup)

      # extract features for Fisher's test
      # a = gene in cluster and in marker database
      a = sum(subVar$gene %in% subMark$gene)
      # b = gene in cluster but not in marker database
      b = length(subVar$gene) - a
      # c = gene in marker database but not in cluster
      c = length(subMark$gene) -a
      # d = leftover from universe
      d = length(geneUniverse) - a - b - c

      # if there are no overlaps go to next iteration
      if (a==0){
        next
      }

      # extract the list of overlapping genes
      overlapGenes = paste(subVar$gene[subVar$gene %in% subMark$gene],
                           collapse = ",")

      # perform Fisher's test
      deTable = matrix(c(a, b, c, d),
                       nrow = 2,
                       dimnames = list(marker=c("yes","no"),
                                       cluster=c("yes","no")))
      FT = fisher.test(deTable, alternative = "greater")

      # put result into a data.frame
      if (systemLevel == "tissue"){

        if (byDatabase){
          tissue = unlist(strsplit(helpGroup, split = "__"))[1]
          database = unlist(strsplit(helpGroup, split = "__"))[2]

          helpTable = data.frame(cluster = clust,
                                 tissue = tissue,
                                 database = database,
                                 p = unname(FT$p.value),
                                 oddsRatio = unname(FT$estimate),
                                 n_olaps = a,
                                 overlaps = overlapGenes)
        } else {
          helpTable = data.frame(cluster = clust,
                                 tissue = helpGroup,
                                 p = unname(FT$p.value),
                                 oddsRatio = unname(FT$estimate),
                                 n_olaps = a,
                                 overlaps = overlapGenes)
        }

      } else if (systemLevel == "cellType"){
        if (byDatabase){
          cellType = unlist(strsplit(helpGroup, split = "__"))[1]
          database = unlist(strsplit(helpGroup, split = "__"))[2]

          helpTable = data.frame(cluster = clust,
                                 cellType = cellType,
                                 database = database,
                                 p = unname(FT$p.value),
                                 oddsRatio = unname(FT$estimate),
                                 n_olaps = a,
                                 overlaps = overlapGenes)
        } else {
          helpTable = data.frame(cluster = clust,
                                 cellType = helpGroup,
                                 p = unname(FT$p.value),
                                 oddsRatio = unname(FT$estimate),
                                 n_olaps = a,
                                 overlaps = overlapGenes)
        }
      }

      # check if variable FisherTable was created (in the first iteration with
      # found overlaps) - if not, put results into it / if yes, just attach
      # new results
      if (!exists(x = "FisherTable")){
        FisherTable = helpTable
      } else {
        FisherTable  = rbind(FisherTable, helpTable)
      }
    }

  }

  # calculate adjusted p-values, sort table by cluster and adjusted p-value
  # calculate -log10(padj)
  FisherTable$padj =  p.adjust(FisherTable$p, method = "BH")
  FisherTable= FisherTable %>%
    arrange(cluster, padj) %>%
    mutate(neg_log10_padj = -log10(padj))

  return(FisherTable)
}

#' #' The function plotMarkerEnrich plots results from \code{markerEnrich} function
#' #
#' #' @param enrichRes An output from \code{markerEnrich} function in form of
#' #'  data.frame
#' #' @return A ggplot object
#' #'
#' #' @export
#' #' @examples
#' #' varMarkerFile = system.file("extdata", "cluster_variable_features.csv",
#' #'                             package = "BoneCellType")
#' #' varMarkers = read.csv(varMarkerFile)
#' #' markerDatabase = markerDatabase
#' #' enrichRes = markerEnrich(varMarkers, markerDatabase)
#' #' plotMarkerEnrich(enrichRes)
#' plotMarkerEnrich = function(enrichRes){
#'
#' }

