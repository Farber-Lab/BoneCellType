#’ Bone marker database
#'
#' Dataset containing bone-related markers.
#'
#' @format A data.frame with 133454 rows and 4 variables
#' \describe{
#'   \item{tissue}{tissue of origin}
#'   \item{cellType}{cell type or origin}
#'   \item{gene}{marker}
#'   \item{database}{source database from which makers were obtained}
#' }
#' @name markerDatabase
#' @source <https://github.com/kkupkova/Mouse-bone-markers/blob/main/BONE_DATABASE.tsv>
"markerDatabase"


#’ Cell type marker database separted by source
#'
#' Dataset (list) containing markers for individual bone-related cell types
#' separated by source database. List names are in format source__cell type.
#' Each item contains associated markers.
#'
#' @format A list with 317 elements
#' @name cellTypeSource_list
#' @source <https://github.com/kkupkova/Mouse-bone-markers/blob/main/BONE_DATABASE_lists/bone_cellType_and_sourceDB.Rdata>
"cellTypeSource_list"


#’ Cell type marker database merged across sources
#'
#' Dataset (list) containing markers for individual bone-related cell types.
#' List names indicate a given cell type. Each item contains associated markers.
#'
#' @format A list with 202 elements
#' @name cellType_list
#' @source <https://github.com/kkupkova/Mouse-bone-markers/blob/main/BONE_DATABASE_lists/bone_cellType.Rdata>
"cellType_list"


#’ Tissue marker database separted by source
#'
#' Dataset (list) containing markers for individual bone-related tissues
#' separated by source database. List names are in format source__tissue.
#' Each item contains associated markers.
#'
#' @format A list with 11 elements
#' @name tissueSource_list
#' @source <https://github.com/kkupkova/Mouse-bone-markers/blob/main/BONE_DATABASE_lists/bone_tissue_and_sourceDB.Rdata>
"tissueSource_list"


#’ Tissue marker database merged across sources
#'
#' Dataset (list) containing markers for individual bone-related tissues.
#' List names indicate a given tissue. Each item contains associated markers.
#'
#' @format A list with 4 elements
#' @name tissue_list
#' @source <https://github.com/kkupkova/Mouse-bone-markers/blob/main/BONE_DATABASE_lists/bone_tissue.Rdata>
"tissue_list"

