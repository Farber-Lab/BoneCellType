# Supporting imports and settings for the package

# You can either use 'import X' or 'importFrom X abcdefg'. importFrom  is
# better practice, but for dplyr and ggplot2 we were simply importing so many functions
# that it makes  sense to just import the whole package
#' @import dplyr
#' @import ggplot2
#' @importFrom stats fisher.test p.adjust
#' @importFrom utils globalVariables
#' @importFrom fgsea fgsea


# Because of some issues with NOTEs on R CMD check and CRAN submission,
# (see here: http://stackoverflow.com/questions/9439256/)
# I have to register stuff used in data.table as non-standard evaluation,
# in order to pass some R CMD check NOTES.
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "avg_log2FC", "cluster", "gene", "myGroup", "padj", "cellTypeSource_list",
    "cellType_list", "tissueSource_list", "tissue_list", "markerDatabase"))
}
