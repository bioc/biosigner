#### biosign (MultiDataSet) ####

#' @rdname biosign
#' @export
setMethod("biosign", signature(x = "MultiDataSet"),
          function(x,
                   y,
                   seedI = NULL,
                   fig.pdfC = c("none", "interactive", "myfile.pdf")[2],                   
                   info.txtC = c("none", "interactive", "myfile.txt")[2],
                   ...) {
            
            if (!(info.txtC %in% c("none", "interactive")))
              sink(info.txtC, append = TRUE)
            
            infTxtC <- info.txtC
            if (infTxtC != "none")
              infTxtC <- "interactive"
            
            if (!(fig.pdfC %in% c("none", "interactive")))
              grDevices::pdf(fig.pdfC)
            
            figPdfC <- fig.pdfC
            if (figPdfC != "none")
              figPdfC <- "none"
            
            biosignMsetLs <- vector(mode = "list",
                                    length = length(names(x)))
            names(biosignMsetLs) <- names(x)
            
            for (setC in names(x)) {
              
              if (info.txtC != "none")
                cat("\n\nSelecting the features for the '", setC, "' dataset:\n", sep = "")
              
              plotL <- TRUE
              
              setBiosign <- tryCatch(biosigner::biosign(x[[setC]],
                                                        y = y,
                                                        fig.pdfC = figPdfC,
                                                        info.txtC = infTxtC,
                                                        ...),
                                     error = function(e) NULL)
              
              if (is.null(setBiosign)) {
                
                setBiosign <- new("biosign")
                setBiosign@eset <- x[[setC]]
                plotL <- FALSE
                
              } else {
                
                tierMC <- setBiosign@tierMC
                
                if (is.null(tierMC) ||
                    sum(apply(setBiosign@tierMC, 2,
                              function(colVc) sum(colVc == "S"))) < 1)
                  plotL <- FALSE
                
              }
              
              biosignMsetLs[[setC]] <- setBiosign
              
              if (fig.pdfC != "none" && plotL) {
                
                plot(setBiosign,
                     typeC = "tier",
                     plotSubC = paste0("[", setC, "]"),
                     ...)
                
                plot(setBiosign,
                     typeC = "boxplot",
                     plotSubC = paste0("[", setC, "]"),
                     ...)
                
              }
              
            }
            
            if (!(fig.pdfC %in% c("none", "interactive")))
              grDevices::dev.off()
            
            if (!(info.txtC %in% c("none", "interactive")))
              sink()
            
            biosignMset <- new("biosignMultiDataSet")
            biosignMset@biosignLs <- biosignMsetLs
            
            return(invisible(biosignMset))
            
          })


#### biosign (ExpressionSet) ####

#' @rdname biosign
#' @export
setMethod("biosign", signature(x = "ExpressionSet"),
          function(x, y, ...) {
            
            if (!(y %in% colnames(Biobase::pData(x)))) {
              stop("'y' must be the name of a column of the phenoData slot of the 'ExpressionSet' object", call. = FALSE)
            } else {
              rspFcVc <- Biobase::pData(x)[, y]
              bsg <- biosign(t(Biobase::exprs(x)), rspFcVc, ...)
            }
            
            tierMC <- bsg@tierMC
            
            fdaDF <- Biobase::fData(x)
            
            if (!is.null(tierMC)) {
              
              for (colI in 1:ncol(tierMC)) {
                
                tierVc <- ropls:::.genVec(x, "feature", "character")
                
                tierVc[rownames(tierMC)] <- tierMC[, colI]
                
                fdaDF[, paste(y,
                              "biosign",
                              gsub("randomforest", "forest", colnames(tierMC)[colI]),
                              sep = "_")] <- tierVc
                
              }
              
            }
            
            Biobase::fData(x) <- fdaDF
            
            bsg@eset <- x
            
            return(invisible(bsg))
            
          })


#### biosign (data.frame) ####

#' @rdname biosign
#' @export
setMethod("biosign", signature(x = "data.frame"),
          function(x, y, ...) {
            if (!all(sapply(x, data.class) == "numeric")) {
              stop("'x' data frame must contain columns of 'numeric' vectors only", call. = FALSE)
            } else
              x <- as.matrix(x)
            bsg <- biosign(x, y, ...)
            
            return(invisible(bsg))
            
          })


####    biosign (matrix)   ####

#' Builds the molecular signature.
#'
#' Main function of the 'biosigner' package. For each of the available
#' classifiers (PLS-DA, Random Forest, and SVM), the significant features are
#' selected and the corresponding models are built.
#'
#' @name biosign
#' @rdname biosign
#' @aliases biosign biosign,data.frame-method biosign,matrix-method
#' biosign,ExpressionSet-method
#' @docType methods
#' @param x Numerical data frame or matrix (observations x variables), or
#' ExpressionSet object with non empty assayData and phenoData; NAs are allowed
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
#' @param printL Logical: deprecated: use the 'info.txtC' argument instead
#' @param plotL Logical: deprecated: use the 'fig.pdfC' argument instead
#' @param .sinkC Character: deprecated: use the 'info.txtC' argument instead
#' @param ... Currently not used.
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
#' @export
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
#' ## Application to an ExpressionSet
#' 
#' diaSet <- ExpressionSet(assayData = t(dataMatrix), 
#'                         phenoData = new("AnnotatedDataFrame", 
#'                                         data = sampleMetadata), 
#'                         featureData = new("AnnotatedDataFrame", 
#'                                           data = variableMetadata),
#'                         experimentData = new("MIAME", 
#'                                              title = "diaplasma"))
#'                                              
#' diaSign <- biosign(diaSet, "type", bootI = 5)
#' diaSet <- getEset(diaSign)
#' head(fData(diaSet))
#' 
#' detach(diaplasma)
#' 
#' ## Application to a MultiDataSet
#' 
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
#' # Building PLS-DA models for the cancer type, and getting back the updated MultiDataSet
#' nciPlsda <- ropls::opls(nciMset, "cancer", predI = 2)
#' nciMset <- ropls::getMset(nciPlsda)
#' # Selecting the significant features for PLS-DA, RF, and SVM classifiers, and getting back the updated MultiDataSet
#' nciBiosign <- biosigner::biosign(nciMset, "cancer")
#' nciMset <- biosigner::getMset(nciBiosign)
setMethod("biosign", signature(x = "matrix"),
          function(x,
                   y,
                   methodVc = c("all", "plsda", "randomforest", "svm")[1],
                   bootI = 50,
                   pvalN = 0.05,
                   
                   permI = 1,
                   fixRankL = FALSE,
                   
                   seedI = 123,
                   plotSubC = NA,
                   fig.pdfC = c("none", "interactive", "myfile.pdf")[2],                   
                   info.txtC = c("none", "interactive", "myfile.txt")[2],                   
                   
                   printL = TRUE,
                   plotL = TRUE,
                   .sinkC = NULL,
                   ...) {
            
            if (!printL) {
              warning("'printL' argument is deprecated; use 'info.txtC' instead.",
                      call. = FALSE)
              info.txtC <- "none"
            }
            
            if (!plotL) {
              warning("'plotL' argument is deprecated; use 'fig.pdfC' instead.",
                      call. = FALSE)
              fig.pdfC <- "none"
            }
            
            if (!is.null(.sinkC)) {
              warning("'.sinkC' argument is deprecated; use 'info.txtC' instead.",
                      call. = FALSE)
              info.txtC <- .sinkC
            }
            
            if (is.null(info.txtC)) {
              warning("'info.txtC = NULL' argument value is deprecated; use 'info.txtC = 'none'' instead.",
                      call. = FALSE)
              info.txtC <- 'none'
            }
            
            if (is.na(info.txtC)) {
              warning("'info.txtC = NA' argument value is deprecated; use 'info.txtC = 'interactive'' instead.",
                      call. = FALSE)
              info.txtC <- 'interactive'
            }
            
            if (is.null(fig.pdfC)) {
              warning("'fig.pdfC = NULL' argument value is deprecated; use 'fig.pdfC = 'none'' instead.",
                      call. = FALSE)
              fig.pdfC <- 'none'
            }
            
            if (is.na(fig.pdfC)) {
              warning("'fig.pdfC = NA' argument value is deprecated; use 'fig.pdfC = 'interactive'' instead.",
                      call. = FALSE)
              fig.pdfC <- 'interactive'
            }
            
            if (!(info.txtC %in% c("none", "interactive")))
              sink(info.txtC, append = TRUE)
            
            if (mode(x) != "numeric")
              stop("'x' matrix must be of 'numeric' mode", call. = FALSE)
            xMN <- x
            
            if (any(apply(xMN, 2, function(colVn) all(is.na(colVn)))))
              stop("'x' contains columns with 'NA' only", call. = FALSE)
            
            if (is.null(colnames(xMN))) {
              colnames(xMN) <- paste0("V", 1:ncol(xMN))
              warning(paste0("Missing column names in 'x' (i.e., variable names): Names are being set to '",
                             head(colnames(xMN), 1),
                             "' ... '",
                             tail(colnames(xMN), 1),
                             "'"),
                      call. = FALSE)
            }
            
            if (is.character(y)) {
              y <- factor(y)
              # warning("'y' character vector converted to a factor with levels '",
              #         paste(levels(y), collapse = "', '"),
              #         "'",
              #         call. = FALSE)
            } else if (!is.factor(y)) {
              stop("'y' must be a character vector or a factor", call. = FALSE)
            }
            
            
            if (length(y) != nrow(xMN)) {
              stop("'y' factor length must be equal to the number of rows of 'x'", call. = FALSE)
            } else if (length(levels(y)) != 2) {
              stop("'y' must have two levels", call. = FALSE)
            } else if (any(is.na(y))) {
              stop("'y' must not contain missing ('NA') values", call. = FALSE)
            } else
              yFc <- y
            
            methodAllVc <- c("plsda", "randomforest", "svm")
            if (length(methodVc) == 1 && methodVc == "all") {
              methodVc <- methodAllVc
            } else if (!all(methodVc %in% methodAllVc))
              stop("The following method(s) is/are not available: '",
                   paste(methodVc[!(methodVc %in% methodAllVc)], collapse = "', '"),
                   "'",
                   call. = FALSE)
            
            if (pvalN < 0 || pvalN > 1)
              stop("The p-value threshold 'pvalN' (",
                   pvalN,
                   ") must be between 0 and 1",
                   call. = FALSE)
            
            
            ## signatures
            
            datasetLs <- list(dataMatrix = xMN,
                              sampleMetadata = data.frame(row.names = rownames(xMN),
                                                          respC = yFc),
                              variableMetadata = data.frame(row.names = colnames(xMN),
                                                            name. = colnames(xMN)))
            tierMN <- matrix(0, nrow = ncol(xMN), ncol = length(methodVc))
            rownames(tierMN) <- colnames(xMN)
            colnames(tierMN) <- methodVc
            
            accuVn <- numeric(length(methodVc))
            stopIterVl <- logical(length(methodVc))
            names(stopIterVl) <- names(accuVn) <- methodVc
            
            for (methodC in methodVc) {
              
              if (!is.null(seedI))
                set.seed(seedI)
              
              optWrnN <- options()$warn
              options(warn = -1)
              
              tierMetLs <- getTierF(datasetLs = datasetLs,
                                    methodC = methodC,
                                    bootI = bootI,
                                    permI = permI,
                                    pvalN = pvalN,
                                    fixRankL = fixRankL)
              
              options(warn = optWrnN)
              
              if (!is.null(seedI))
                set.seed(NULL)
              
              tierMN[, methodC] <- tierMetLs[["tierVn"]]
              accuVn[methodC] <- tierMetLs[["accuracyN"]]
              stopIterVl[methodC] <- tierMetLs[["stopIterL"]]
              
            }
            
            ## Formatting, ordering, and displaying results
            
            xSubMNAS <- xSubMN <- signatureLsAS <- signatureLs <- modelLsAS <- modelLs <- tierMC <- NULL
            
            ## Accuracy of full, S+A and S models
            
            accuracyMN <- matrix(NA, nrow = 3, ncol = length(methodVc),
                                 dimnames = list(c("Full", "AS", "S"), methodVc))
            
            modelLs <- vector(mode = "list", length = length(methodVc))
            names(modelLs) <- methodVc
            signatureLs <- modelLsAS <- modelLs ## void models are empty lists
            for (sgnI in 1:length(signatureLs))
              signatureLs[[sgnI]] <- character() ## void signatures are character(0)
            signatureLsAS <- signatureLs
            
            ## translating tiers (building tierMC)
            
            if (max(tierMN) > 0) {
              for (methodC in methodVc) {
                tierVn <- tierMN[, methodC]
                translate <- rep("E", length(unique(tierVn)))
                if (max(tierVn) > 0) {
                  ## if(stop.rec[[methodC]]) {
                  if (stopIterVl[methodC]) {
                    translate[max(tierVn) + 1] <- "S"
                    if (max(tierVn) > 1)
                      translate[max(tierVn)] <- "A"
                    if (max(tierVn) > 2)
                      translate[max(tierVn) - 1] <- "B"
                    if (max(tierVn) > 3)
                      translate[max(tierVn) - 2] <- "C"
                    if (max(tierVn) > 4)
                      translate[max(tierVn) - 3] <- "D"
                  }
                  else{
                    translate[max(tierVn) + 1] <- "A"
                    if (max(tierVn) > 1)
                      translate[max(tierVn)] <- "B"
                    if (max(tierVn) > 2)
                      translate[max(tierVn) - 1] <- "C"
                    if (max(tierVn) > 3)
                      translate[max(tierVn) - 2] <- "D"
                  }
                }
                tierVn <- translate[tierVn + 1]
                tierMC <- cbind(tierMC, tierVn)
              }
              dimnames(tierMC) <- dimnames(tierMN)
              
              ## ordering tiers (ordering tierMC)
              
              tierFullVc <- c("S", LETTERS[1:5])
              tierFullVn <- length(tierFullVc):1
              names(tierFullVn) <- tierFullVc
              tierOrdMC <- NULL
              for (tierC in tierFullVc) {
                if (any(tierC %in% c(tierMC))) {
                  
                  rowSelVl <- rowSums(tierMC == tierC) > 0
                  tierMaxMC <- tierMC[rowSelVl, , drop = FALSE]
                  tierMaxMN <- tierMaxMC
                  for (j in 1:ncol(tierMaxMC))
                    tierMaxMN[, j] <- tierFullVn[tierMaxMC[, j]]
                  mode(tierMaxMN) <- "numeric"
                  tierMaxMC <- tierMaxMC[order(rowSums(tierMaxMN), decreasing = TRUE), , drop = FALSE]
                  tierOrdMC <- rbind(tierOrdMC, tierMaxMC)
                  tierMC <- tierMC[!rowSelVl, , drop = FALSE]
                }
              }
              tierMC <- tierOrdMC
              
              rm(tierMN)
              
              ## Accuracy of full, S+A and S models
              
              for (methodC in methodVc) {
                ## accuracyMN["Full", methodC] <- fsiResLs[[methodC]]$accuracyN
                accuracyMN["Full", methodC] <- accuVn[methodC]
                
                ## 'S+A' outputs
                datasetMethodLs <- datasetLs
                
                signatureLsAS[[methodC]] <- rownames(tierMC)[which(tierMC[, methodC] %in% c("S","A"))]
                
                if (length(signatureLsAS[[methodC]])) {
                  datasetMethodLs$variableMetadata <- datasetMethodLs$variableMetadata[signatureLsAS[[methodC]], , drop = FALSE]
                  datasetMethodLs$dataMatrix  <- datasetMethodLs$dataMatrix[ , signatureLsAS[[methodC]], drop = FALSE]
                  
                  if (!is.null(seedI))
                    set.seed(seedI)
                  
                  optWrnN <- options()$warn
                  options(warn = -1)

                  my.fsi <- getBootSignificanceF(datasetMethodLs,
                                                 "respC",
                                                 methodC,
                                                 bootI = bootI,
                                                 permI = 0,
                                                 fullModelL = TRUE)
                  
                  options(warn = optWrnN)
                  
                  if (!is.null(seedI))
                    set.seed(NULL)
                  
                  accuracyMN["AS", methodC] <- my.fsi$accuracyN
                  modelLsAS[[methodC]] <- my.fsi$model
                }
                
                ## 'S' outputs
                datasetMethodLs <- datasetLs
                
                signatureLs[[methodC]] <- rownames(tierMC)[which(tierMC[, methodC] == "S")]
                
                if (length(signatureLs[[methodC]])) {
                  datasetMethodLs$variableMetadata <- datasetMethodLs$variableMetadata[signatureLs[[methodC]], , drop = FALSE]
                  datasetMethodLs$dataMatrix  <- datasetMethodLs$dataMatrix[ , signatureLs[[methodC]], drop = FALSE]
                  
                  if (!is.null(seedI))
                    set.seed(seedI)
                  
                  optWrnN <- options()$warn
                  options(warn = -1)
                  
                  my.fsi <- getBootSignificanceF(datasetMethodLs,
                                                 "respC",
                                                 methodC,
                                                 bootI = bootI,
                                                 permI = 0,
                                                 fullModelL = TRUE)
                  
                  options(warn = optWrnN)
                  
                  if (!is.null(seedI))
                    set.seed(NULL)
                  
                  accuracyMN["S", methodC] <- my.fsi$accuracyN
                  modelLs[[methodC]] <- my.fsi$model
                }
              }
              signatureLs[["complete"]] <- rownames(tierMC)[apply(tierMC, 1, function(rowVc) sum(rowVc == "S") > 0)]
              xSubMN <- xMN[, signatureLs[["complete"]], drop = FALSE] ## for boxplotting
              signatureLsAS[["complete"]] <- rownames(tierMC)[apply(tierMC, 1, function(rowVc) sum(rowVc %in% c("S", "A")) > 0)]
              xSubMNAS <- xMN[, signatureLsAS[["complete"]], drop = FALSE] ## for boxplotting
              
            } else {
              
              tierMC <- matrix("E", nrow = nrow(tierMN), ncol = ncol(tierMN),
                               dimnames = dimnames(tierMN))
              
              ## Accuracy of full model
              
              for (methodC in methodVc) {
                
                accuracyMN["Full", methodC] <- accuVn[methodC]
                
                xSubMN <- xMN[, signatureLs[["complete"]], drop = FALSE] ## for boxplotting
                
                xSubMNAS <- xMN[, signatureLsAS[["complete"]], drop = FALSE] ## for boxplotting
                
              }
              
            }
            
            AS <- list(modelLs = modelLsAS,
                       signatureLs = signatureLsAS,
                       xSubMN = xSubMNAS)
            
            bsg <- new("biosign")
            bsg@methodVc <- methodVc
            bsg@accuracyMN <- accuracyMN
            bsg@tierMC <- tierMC
            bsg@yFc <- yFc
            bsg@modelLs <- modelLs
            bsg@signatureLs <- signatureLs
            bsg@xSubMN <- xSubMN
            bsg@AS <- AS
            
            ## Printing
            
            if (info.txtC != "none") {
              show(bsg)
            }
            
            ## Plotting
            
            if (!all(is.na(bsg@accuracyMN["S", ])) && fig.pdfC != "none")
              plot(bsg, typeC = "tier", plotSubC = plotSubC, fig.pdfC = fig.pdfC)
            
            ## Closing connection
            
            if (!(info.txtC %in% c("none", "interactive")))
              sink()
            
            ## Returning
            
            return(invisible(bsg))
            
            
          })

####    predict    ####

#' Predict method for 'biosign' signature objects
#'
#' This function predicts values based upon \code{biosign} classifiers trained
#' on the 'S' signature
#'
#'
#' @aliases predict.biosign predict,biosign-method
#' @param object An S4 object of class \code{biosign}, created by
#' \code{biosign} function.
#' @param newdata Either a data frame or a matrix, containing numeric columns
#' only, with column names identical to the 'x' used for model training with
#' 'biosign'.
#' @param tierMaxC Character: Maximum level of tiers to display: Either 'S'or
#' 'A'.
#' @param ... Currently not used.
#' @return Data frame with the predictions for each classifier as factor
#' columns.
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
#' ## training the classifiers
#' ## a bootI = 5 number of bootstraps is used for this example
#' ## we recommend to keep the default bootI = 50 value for your analyzes
#'
#' set.seed(123)
#' diaSign <- biosign(dataMatrix, sampleMetadata[, "type"], bootI = 5)
#'
#' ## fitted values (for the subsets restricted to the 'S' signatures)
#' sFitDF <- predict(diaSign)
#'
#' ## confusion tables
#' print(lapply(sFitDF, function(predFc) table(actual = sampleMetadata[,
#' "type"], predicted = predFc)))
#'
#' ## balanced accuracies
#' sapply(sFitDF, function(predFc) { conf <- table(sampleMetadata[,
#' "type"], predFc)
#' conf <- sweep(conf, 1, rowSums(conf), "/")
#' mean(diag(conf))
#' })
#' ## note that these values are slightly different from the accuracies
#' ## returned by biosign because the latter are computed by using the
#' ## resampling scheme selected by the bootI or crossvalI arguments
#' getAccuracyMN(diaSign)["S", ]
#'
#' detach(diaplasma)
#'
#' @rdname predict
#' @export
setMethod("predict", signature(object = "biosign"),
          function(object, newdata, tierMaxC = "S", ...)
          {
            
            if (any(!(tierMaxC %in% c("S", "A"))))
              stop("'tierMaxC' argument must be 'S' or 'A' for predictions", call. = FALSE)
            
            switch(tierMaxC,
                   S = {
                     signatureLs <- object@signatureLs
                     modelLs <- object@modelLs
                   },
                   A = {
                     signatureLs <- object@AS[["signatureLs"]]
                     modelLs <- object@AS[["modelLs"]]
                   })
            
            if (missing(newdata)) {
              
              fitLs <- lapply(modelLs,
                              function(model) {
                                if(is.null(model))
                                  return(NULL)
                                else
                                  return(switch(class(model),
                                                opls = fitted(model),
                                                randomForest = model[["predicted"]],
                                                svm = model[["fitted"]]))
                              })
              
              fitLs <- fitLs[!sapply(fitLs, is.null)]
              
              if (length(fitLs) == 0) {
                warning("Empty signatures for all classifiers up to tier '", tierMaxC, "'; fitted output is set to NULL", call. = FALSE)
                return(invisible(NULL))
              } else
                return(as.data.frame(fitLs))
              
            } else {
              
              if (is.data.frame(newdata)) {
                if(!all(sapply(newdata, data.class) == "numeric")) {
                  stop("'newdata' data frame must contain numeric columns only", call. = FALSE)
                }
              } else if(is.matrix(newdata)) {
                if(mode(newdata) != "numeric") {
                  stop("'newdata' matrix must be of 'numeric' mode", call. = FALSE)
                } else
                  newdata <- as.data.frame(newdata)
              } else
                stop("'newdata' must be either a data.frame or a matrix", call. = FALSE)
              
              if(is.null(colnames(newdata))) {
                colnames(newdata) <- paste0("V", 1:ncol(newdata))
                warning(paste0("Missing column names in 'newdata' (i.e., variable names): Names are being set to '",
                               head(colnames(newdata), 1),
                               "' ... '",
                               tail(colnames(newdata), 1),
                               "'"),
                        call. = FALSE)
              }
              
              if (length(signatureLs[["complete"]]) == 0){
                
                warning("Signatures from all classifiers up to tier '", tierMaxC, "' are empty; prediction output is set to NULL", call. = FALSE)
                
                return(invisible(NULL))
                
              } else {
                
                predDF <- as.data.frame(lapply(object@methodVc,
                                               function(methodC) {
                                                 signVc <- signatureLs[[methodC]]
                                                 if(length(signVc)) {
                                                   signOutVc <- signVc[!(signVc %in% colnames(newdata))]
                                                   if(length(signOutVc)) {
                                                     stop("The following variables from the ", methodC, " ", switch(tierMaxC, S = "S", A = "S+A"), " signature are not found as column names in 'newdata': '", paste(signOutVc, collapse = "', '"), "'.", call. = FALSE)
                                                   } else
                                                     return(predict(modelLs[[methodC]],
                                                                    newdata = newdata[, signVc, drop = FALSE]))
                                                 } else
                                                   return(rep(NA, nrow(newdata))) ## /!\ logical mode
                                               }))
                colnames(predDF) <- object@methodVc
                predDF <- predDF[, !apply(predDF, 2, function(predVn) all(is.na(predVn))), drop = FALSE]
                
                return(predDF)
                
              }
              
            }
            
          })


####    show    ####

#' Show method for 'biosign' signature objects
#'
#' Prints the selected features and the accuracies of the classifiers.
#'
#' @aliases show.biosign show,biosign-method
#' @param object An S4 object of class \code{biosign}, created by the
#' \code{biosign} function.
#' @return Invisible.
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
#' diaSign
#'
#' detach(diaplasma)
#'
#' @rdname show
#' @export
setMethod("show",
          "biosign",
          function(object)
          {
            
            tierMaxC <- "S"
            
            ## if(is.null(object@accuracyMN)) {
            if(all(is.na(object@accuracyMN["S", ]))) {
              
              cat("No significant variable found for the selected classifier(s): '", paste(object@methodVc, collapse = "', '"), "'\n", sep = "")
              
            } else {
              
              tierFullVc <- c("S", "A", "B", "C", "D", "E")
              
              if (any(!(tierMaxC %in% tierFullVc))) {
                stop("'tierMaxC' argument must be in '", paste(tierFullVc, collapse = "', '"), "'", call. = FALSE)
              } else
                tierVc <- tierFullVc[1:which(tierFullVc == tierMaxC)]
              
              if (sum(object@tierMC %in% tierVc)) {
                if (length(setdiff(tierVc, c(object@tierMC)))) {
                  tierNotFoundVc <- tierFullVc[tierFullVc %in% setdiff(tierVc, c(object@tierMC))]
                  warning("tierMC does not contain the following values: '", paste(tierNotFoundVc, collapse = "', '"), "'", call. = FALSE)
                }
                
                cat("Significant features from '", paste(tierVc, collapse = "', '"), "' groups:\n", sep = "")
                print(object@tierMC[apply(object@tierMC, 1, function(rowVc) sum(rowVc %in% tierVc) > 0), ,
                                    drop = FALSE])
                
                cat("Accuracy:\n", sep = "")
                print(round(object@accuracyMN, 3))
                
                invisible(NULL)
                
              }
              
            }
            
          })


####    getEset    ####

#' getEset method
#'
#' Extracts the complemented ExpressionSet when biosign has been applied to an ExpressionSet
#'
#' @aliases getEset getEset, biosign-method
#' @param object An S4 object of class \code{biosign}, created by \code{biosign}
#' function.
#' @return An S4 object of class \code{ExpressionSet} which contains the dataMatrix (t(exprs(eset))),
#' and the sampleMetadata (pData(eset)) and variableMetadata (fData(eset)) with the additional columns
#' containing the computed tiers for each feature and each classifier.
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#' 
#' ## loading the diaplasma dataset
#'
#' data(diaplasma)
#' attach(diaplasma)
#' 
#' ## building the ExpresssionSet instance
#' 
#' diaSet <- Biobase::ExpressionSet(assayData = t(dataMatrix), 
#'                                  phenoData = new("AnnotatedDataFrame", 
#'                                                  data = sampleMetadata), 
#'                                  featureData = new("AnnotatedDataFrame", 
#'                                                    data = variableMetadata),
#'                                  experimentData = new("MIAME", 
#'                                                       title = "diaplasma"))
#'
#' ## restricting to a smaller dataset for this example
#'
#' featureSelVl <- variableMetadata[, "mzmed"] >= 490 & variableMetadata[, "mzmed"] < 500
#' diaSet <- diaSet[featureSelVl, ]
#'
#' ## signature selection for all 3 classifiers
#' ## a bootI = 5 number of bootstraps is used for this example
#' ## we recommend to keep the default bootI = 50 value for your analyzes
#'
#' set.seed(123)
#' diaSign <- biosign(diaSet, "type", bootI = 5)
#' 
#' diaSet <- biosigner::getEset(diaSign)
#' head(Biobase::pData(diaSet))
#' head(Biobase::fData(diaSet))
#' 
#' detach(diaplasma)
#'
#' @rdname getEset
#' @export
setMethod("getEset", "biosign",
          function(object) {
            return(object@eset)
          })


####     getAccuracyMN    ####

#' Accuracies of the full model and the models restricted to the signatures
#'
#' Balanced accuracies for the full models, and the models restricted to the
#' 'S' and 'AS' signatures
#'
#' @aliases getAccuracyMN getAccuracyMN,biosign-method
#' @param object An S4 object of class \code{biosign}, created by the
#' \code{biosign} function.
#' @param ... Currently not used.
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
setMethod("getAccuracyMN", "biosign",
          function(object) {
            return(object@accuracyMN)
          })


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
#' @param ... Currently not used.
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
setMethod("getSignatureLs", "biosign",
          function(object, tierC = c("S", "AS")[1]) {
            if(tierC == "S") {
              return(object@signatureLs)
            } else if(tierC == "AS") {
              return(object@AS$signatureLs)
            } else
              stop("'tierC' argument must be either 'S' or 'AS'")
          })




