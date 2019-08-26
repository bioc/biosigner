####    getMset (biosignMultiDataSet)   ####

#' getMset method
#'
#' Extracts the complemented MultiDataSet when biosign has been applied to a MultiDataSet
#'
#' @aliases getMset getMset, biosignMultiDataSet-method
#' @param object An S4 object of class \code{biosignMultiDataSet}, created by \code{biosign}
#' function applied to a MultiDataSet
#' @param ... Currently not used
#' @return An S4 object of class \code{MultiDataSet}.
#' @rdname getMset
#' @export
#' @examples
#' # Loading the 'NCI60_4arrays' from the 'omicade4' package
#' data("NCI60_4arrays", package = "omicade4")
#' # Selecting two of the four datasets
#' setNamesVc <- c("agilent", "hgu95")
#' # Creating the MultiDataSet instance
#' nciMset <- MultiDataSet::createMultiDataSet()
#' # Adding the two datasets as ExpressionSet instances
#' for (setC in setNamesVc) {
#'   # Getting the data
#'   exprMN <- as.matrix(NCI60_4arrays[[setC]])
#'   pdataDF <- data.frame(row.names = colnames(exprMN),
#'                         cancer = substr(colnames(exprMN), 1, 2),
#'                         stringsAsFactors = FALSE)
#'   fdataDF <- data.frame(row.names = rownames(exprMN),
#'                         name = rownames(exprMN),
#'                         stringsAsFactors = FALSE)
#'   # Building the ExpressionSet
#'   eset <- Biobase::ExpressionSet(assayData = exprMN,
#'                                  phenoData = new("AnnotatedDataFrame",
#'                                                  data = pdataDF),
#'                                  featureData = new("AnnotatedDataFrame",
#'                                                    data = fdataDF),
#'                                  experimentData = new("MIAME",
#'                                                       title = setC))
#'   # Adding to the MultiDataSet
#'   nciMset <- MultiDataSet::add_eset(nciMset, eset, dataset.type = setC,
#'                                     GRanges = NA, warnings = FALSE)
#' }
#' # Restricting to the 'ME' and 'LE' cancer types
#' sampleNamesVc <- Biobase::sampleNames(nciMset[["agilent"]])
#' cancerTypeVc <- Biobase::pData(nciMset[["agilent"]])[, "cancer"]
#' nciMset <- nciMset[sampleNamesVc[cancerTypeVc %in% c("ME", "LE")], ]
#' # Summary of the MultiDataSet
#' nciMset
#' # Selecting the significant features for PLS-DA, RF, and SVM classifiers, and getting back the updated MultiDataSet
#' nciBiosign <- biosigner::biosign(nciMset, "cancer")
#' nciMset <- biosigner::getMset(nciBiosign)
#' # In the updated MultiDataSet, the updated featureData now contains the cancer_biosign_'classifier' columns
#' indicating the selected features
#' lapply(fData(nciMset), head)
setMethod("getMset", "biosignMultiDataSet",
          function(object) {
            Mset <- MultiDataSet::createMultiDataSet()
            for (setI in 1:length(object@biosignLs)) {
              Mset <- MultiDataSet::add_eset(Mset,
                                             ropls::getEset(object@biosignLs[[setI]]),
                                             dataset.type = names(object@biosignLs)[setI],
                                             GRanges = NA,
                                             overwrite = TRUE,
                                             warnings = FALSE)
            }
            return(Mset)
          })