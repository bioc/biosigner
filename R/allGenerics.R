#### biosign ####

#' Builds the molecular signature.
#'
#' Main function of the 'biosigner' package. For each of the available
#' classifiers (PLS-DA, Random Forest, and SVM), the significant features are
#' selected and the corresponding models are built.
#'
#' @name biosign
#' @aliases biosign biosign,matrix-method biosign,data.frame-method 
#' biosign,SummarizedExperiment-method biosign,ExpressionSet-method
#' biosign,MultiAssayExperiment-method biosign,MultiDataSet_method
#' @docType methods
#' @param x Numerical data frame or matrix (observations x variables), or
#' SummarizedExperiment (or ExpressionSet) object ; NAs are allowed
#' for PLS-DA but for SVM, samples with NA will be removed
#' @param y Two-level factor corresponding to the class labels, or a character
#' indicating the name of the column of the pData to be used, when x is an
#' ExpressionSet object
#' @param methodVc Character vector: Either one or all of the following
#' classifiers: Partial Least Squares Discriminant Analysis ('plsda'), or
#' Random Forest ('randomforest'), or Support Vector Machine ('svm')
#' @param bootI Integer: Number of bootstaps for resampling
#' @param pvalN Numeric: To speed up the selection, only variables which
#' significantly improve the model up to two times this threshold (to take into
#' account potential fluctuations) are computed
#' @param permI Integer: Random permutation are used to assess the significance
#' of each new variable included into the model (forward selection)
#' @param fixRankL Logical: Should the initial ranking be computed with the
#' full model only, or as the median of the ranks from the models built on the
#' sampled dataset?
#' @param seedI integer: optional seed to obtain exactly the same signature when
#' rerunning biosigner; default is '123'; set to NULL to prevent seed setting
#' @param plotSubC Character: Graphic subtitle
#' @param fig.pdfC Character: File name with '.pdf' extension for the figure;
#' if 'interactive' (default), figures will be displayed interactively; if 'none',
#' no figure will be generated
#' @param info.txtC Character: File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return An S4 object of class 'biosign' containing the following slots: 1)
#' 'methodVc' character vector: selected classifier(s) ('plsda',
#' 'randomforest', and/or 'svm'), 2) 'accuracyMN' numeric matrix: balanced
#' accuracies for the full models, and the models restricted to the 'S' and
#' 'AS' signatures (predictions are obtained by using the resampling scheme
#' selected with the 'bootI' and 'crossvalI' arguments), 3) 'tierMC' character
#' matrix: contains the tier ('S', 'A', 'B', 'C', 'D', or 'E') of each feature
#' for each classifier (features with tier 'S' have been found significant in
#' all backward selections; features with tier 'A' have been found significant
#' in all but the last selection, and so on), 4) modelLs list: selected
#' classifier(s) trained on the subset restricted to the 'S' features, 5)
#' signatureLs list: 'S' signatures for each classifier; and 6) 'AS' list: 'AS'
#' signatures and corresponding trained classifiers, in addition to the dataset
#' restricted to tiers 'S' and 'A' ('xMN') and the labels ('yFc')
#' @author Philippe Rinaudo and Etienne Thevenot (CEA)
#' @seealso \code{\link{predict.biosign}}, \code{\link{plot.biosign}}
#' @examples
#'
#' ## loading the diaplasma dataset
#'
#' data(diaplasma)
#' attach(diaplasma)
#'
#' ## restricting to a smaller dataset for this example
#'
#' featureSelVl <- variableMetadata[, "mzmed"] >= 490 & variableMetadata[, "mzmed"] < 500
#' dataMatrix <- dataMatrix[, featureSelVl]
#' variableMetadata <- variableMetadata[featureSelVl, ]
#'
#' # signature selection for all 3 classifiers
#' # a bootI = 5 number of bootstraps is used for this example
#' # we recommend to keep the default bootI = 50 value for your analyzes
#'
#' diaSign <- biosign(dataMatrix, sampleMetadata[, "type"], bootI = 5)
#' 
#' ## Application to a SummarizedExperiment
#' 
#' diaplasma.se <- SummarizedExperiment::SummarizedExperiment(assays = list(diaplasma = t(dataMatrix)),
#'                                                            colData = sampleMetadata,
#'                                                            rowData = variableMetadata)
#'                                                            
#' # restricting to the first 100 features to speed up the example
#' 
#' diaplasma.se <- diaplasma.se[1:100, ]
#' 
#' diaplasma.se <- biosign(diaplasma.se, "type", bootI = 5)
#' 
#' head(SummarizedExperiment::rowData(diaplasma.se))
#' 
#' # getting the biosign output
#' 
#' diaplasma_type.biosign <- getBiosign(diaplasma.se)[["type_plsda.forest.svm"]]
#' 
#' getAccuracyMN(diaplasma_type.biosign)
#'
#' ## Application to an ExpressionSet
#' 
#' diaSet <- Biobase::ExpressionSet(assayData = t(dataMatrix), 
#'                                  phenoData = new("AnnotatedDataFrame", 
#'                                                 data = sampleMetadata), 
#'                                  featureData = new("AnnotatedDataFrame", 
#'                                                 data = variableMetadata),
#'                                  experimentData = new("MIAME", 
#'                                                title = "diaplasma"))
#'                                              
#' # restricting to the first 100 features to speed up the example
#' 
#' diaSet <- diaSet[1:100, ]
#'                                              
#' diaSign <- biosign(diaSet, "type", bootI = 5)
#' diaSet <- getEset(diaSign)
#' head(Biobase::fData(diaSet))
#' 
#' detach(diaplasma)
#' 
#' ## Application to a MultiAssayExperiment
#' 
#' data("NCI60", package = "ropls")
#' nci.mae <- NCI60[["mae"]]
#' library(MultiAssayExperiment)
#' 
#' # Cancer types
#' 
#' table(nci.mae$cancer)
#' 
#' # Restricting to the 'ME' and 'LE' cancer types and to the 'agilent' and 'hgu95' datasets
#' 
#' nci.mae <- nci.mae[, nci.mae$cancer %in% c("ME", "LE"), c("agilent", "hgu95")]
#' 
#' # Selecting the significant features for PLS-DA, RF, and SVM classifiers
#' 
#' nci.mae <- biosign(nci.mae, "cancer", bootI = 5)
#' 
#' # Getting the tiers
#' 
#' SummarizedExperiment::rowData(nci.mae[["agilent"]])
#' 
#' # Getting the models
#' 
#' mae_biosign.ls <- getBiosign(nci.mae)
#' 
#' # Name of the models stored in the (metadata of) each SummarizedExperiment object
#' 
#' names(mae_biosign.ls[["agilent"]])
#' 
#' # Visualizing the individual results
#' 
#' for (set.c in names(mae_biosign.ls))
#'   plot(mae_biosign.ls[[set.c]][["cancer_plsda.forest.svm"]],
#'        typeC = "tier",
#'        plotSubC = set.c)
#' 
#' ## Application to a MultiDataSet
#' 
#' data("NCI60", package = "ropls")
#' nci.mds <- NCI60[["mds"]]
#' 
#' # Restricting to the "agilent" and "hgu95" datasets
#' 
#' nci.mds <- nci.mds[, c("agilent", "hgu95")]
#' 
#' # Restricting to the 'ME' and 'LE' cancer types
#' 
#' library(Biobase)
#' sample_names.vc <- Biobase::sampleNames(nci.mds[["agilent"]])
#' cancer_type.vc <- Biobase::pData(nci.mds[["agilent"]])[, "cancer"]
#' nci.mds <- nci.mds[sample_names.vc[cancer_type.vc %in% c("ME", "LE")], ]
#' 
#' # Selecting the significant features for PLS-DA, RF, and SVM classifiers
#' 
#' nci_cancer.biosign <- biosign(nci.mds, "cancer", bootI = 5)
#' 
#' # Getting back the updated MultiDataSet
#' 
#' nci.mds <- getMset(nci_cancer.biosign)
#' 
#' @rdname biosign
#' @export
setGeneric("biosign", function(x,
                               y,
                               methodVc = c("all", "plsda", "randomforest", "svm")[1],
                               bootI = 50,
                               pvalN = 0.05,
                               
                               permI = 1,
                               fixRankL = FALSE,
                               
                               seedI = 123,
                               plotSubC = NA,
                               fig.pdfC = c("none", "interactive", "myfile.pdf")[2],                   
                               info.txtC = c("none", "interactive", "myfile.txt")[2]) standardGeneric("biosign"))


####     getAccuracyMN    ####

#' Accuracies of the full model and the models restricted to the signatures
#'
#' Balanced accuracies for the full models, and the models restricted to the
#' 'S' and 'AS' signatures
#'
#' @aliases getAccuracyMN getAccuracyMN,biosign-method
#' @param object An S4 object of class \code{biosign}, created by the
#' \code{biosign} function.
#' @return A numeric matrix containing the balanced accuracies for the full
#' models, and the models restricted to the 'S' and 'AS' signatures
#' (predictions are obtained by using the resampling scheme selected with the
#' 'bootI' and 'crossvalI' arguments)
#' @author Philippe Rinaudo and Etienne Thevenot (CEA)
#' @examples
#'
#' ## loading the diaplasma dataset
#'
#' data(diaplasma)
#' attach(diaplasma)
#'
#' ## restricting to a smaller dataset for this example
#'
#' featureSelVl <- variableMetadata[, "mzmed"] >= 490 & variableMetadata[, "mzmed"] < 500
#' dataMatrix <- dataMatrix[, featureSelVl]
#' variableMetadata <- variableMetadata[featureSelVl, ]
#'
#' ## signature selection for all 3 classifiers
#' ## a bootI = 5 number of bootstraps is used for this example
#' ## we recommend to keep the default bootI = 50 value for your analyzes
#'
#' set.seed(123)
#' diaSign <- biosign(dataMatrix, sampleMetadata[, "type"], bootI = 5)
#'
#' ## individual boxplot of the selected signatures
#'
#' getAccuracyMN(diaSign)
#' 
#' detach(diaplasma)
#'
#' @rdname getAccuracyMN
#' @export
setGeneric("getAccuracyMN", function(object) {standardGeneric("getAccuracyMN")})


####   getBiosign    ####

#' Getting the biosigner signature from the SummarizedExperiment object
#'
#' The models are extracted as a list
#'
#' @param object An S4 object of class \code{SummarizedExperiment}, once processed by the \code{biosign}
#' method
#' @return List of biosigner outputs contained in the SummarizedExperiment object
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' # Getting the diaplasma data set as a SummarizedExperiment
#' 
#' data(diaplasma)
#' 
#' diaplasma.se <- SummarizedExperiment::SummarizedExperiment(assays = list(diaplasma = t(diaplasma[["dataMatrix"]])),
#'                                                            colData = diaplasma[["sampleMetadata"]],
#'                                                            rowData = diaplasma[["variableMetadata"]])
#' 
#' diaplasma.se <- diaplasma.se[1:100, ]
#' 
#' # Selecting the features
#' 
#' diaplasma.se <- biosign(diaplasma.se, "type", bootI = 5, fig.pdfC = "none")
#' 
#' # Getting the signatures
#' 
#' diaplasma.biosign <- getBiosign(diaplasma.se)[["type_plsda.forest.svm"]]
#' 
#' diaplasma.biosign
#' 
#' @rdname getBiosign
#' @export
setGeneric("getBiosign",
           function(object) {standardGeneric("getBiosign")})


####    getSignatureLs    ####

#' Signatures selected by the models
#'
#' List of 'S' (or 'S' and 'A') signatures for each classifier
#'
#' @aliases getSignatureLs getSignatureLs,biosign-method
#' @param object An S4 object of class \code{biosign}, created by the
#' \code{biosign} function.
#' @param tierC Character: defines whether signatures from the 'S' tier only
#' (default) or the ('S' and 'A') tiers should be returned
#' @return List of 'S' (or 'S' and 'A') signatures for each classifier
#' @author Philippe Rinaudo and Etienne Thevenot (CEA)
#' @examples
#'
#' ## loading the diaplasma dataset
#'
#' data(diaplasma)
#' attach(diaplasma)
#'
#' ## restricting to a smaller dataset for this example
#'
#' featureSelVl <- variableMetadata[, "mzmed"] >= 490 & variableMetadata[, "mzmed"] < 500
#' dataMatrix <- dataMatrix[, featureSelVl]
#' variableMetadata <- variableMetadata[featureSelVl, ]
#'
#' ## signature selection for all 3 classifiers
#' ## a bootI = 5 number of bootstraps is used for this example
#' ## we recommend to keep the default bootI = 50 value for your analyzes
#'
#' set.seed(123)
#' diaSign <- biosign(dataMatrix, sampleMetadata[, "type"], bootI = 5)
#'
#' ## individual boxplot of the selected signatures
#'
#' getSignatureLs(diaSign)
#'
#' detach(diaplasma)
#'
#' @rdname getSignatureLs
#' @export
setGeneric("getSignatureLs", function(object, tierC = c("S", "AS")[1]) {standardGeneric("getSignatureLs")})
