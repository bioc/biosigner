####    biosign (matrix)   ####

#' @rdname biosign
#' @export
setMethod("biosign", signature(x = "matrix"),
          function(x,
                   y,
                   methodVc = c("all", "plsda", "randomforest", "svm")[1],
                   bootI = 50,
                   pvalN = 0.05,
                   
                   permI = 1,
                   fixRankL = FALSE,
                   
                   seedI = 123,
                   plotSubC = "",
                   fig.pdfC = c("none", "interactive", "myfile.pdf")[2],                   
                   info.txtC = c("none", "interactive", "myfile.txt")[2]) {
           
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
              
              if (info.txtC != "none")
                cat("Selecting features for the ", methodC, " model\n", sep = "")
              
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


#### biosign (data.frame) ####

#' @rdname biosign
#' @export
setMethod("biosign", signature(x = "data.frame"),
          function(x,
                   y,
                   methodVc = c("all", "plsda", "randomforest", "svm")[1],
                   bootI = 50,
                   pvalN = 0.05,
                   
                   permI = 1,
                   fixRankL = FALSE,
                   
                   seedI = 123,
                   plotSubC = "",
                   fig.pdfC = c("none", "interactive", "myfile.pdf")[2],                   
                   info.txtC = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!all(sapply(x, data.class) == "numeric")) {
              stop("'x' data frame must contain columns of 'numeric' vectors only", call. = FALSE)
            } else
              x <- as.matrix(x)
            
            bsg <- biosign(x = x,
                           y = y,
                           methodVc = methodVc,
                           bootI = bootI,
                           pvalN = pvalN,
                           
                           permI = permI,
                           fixRankL = fixRankL,
                           
                           seedI = seedI,
                           plotSubC = plotSubC,
                           fig.pdfC = fig.pdfC,                   
                           info.txtC = info.txtC)
            
            return(invisible(bsg))
            
          })


#### biosign (SummarizedExperiment) ####

#' @rdname biosign
#' @export
setMethod("biosign", signature(x = "SummarizedExperiment"),
          function(x,
                   y,
                   methodVc = c("all", "plsda", "randomforest", "svm")[1],
                   bootI = 50,
                   pvalN = 0.05,
                   
                   permI = 1,
                   fixRankL = FALSE,
                   
                   seedI = 123,
                   plotSubC = "",
                   fig.pdfC = c("none", "interactive", "myfile.pdf")[2],                   
                   info.txtC = c("none", "interactive", "myfile.txt")[2]) {
            
            if (length(y) == 1 && is.character(y)) {
              
              if (!(y %in% colnames(SummarizedExperiment::colData(x))))
                stop("'y' must be the name of a column of the colData slot of the 'SummarizedExperiment' object")
              
              bsg <- biosign(x = t(SummarizedExperiment::assay(x)),
                             y = SummarizedExperiment::colData(x)[, y],
                             methodVc = methodVc,
                             bootI = bootI,
                             pvalN = pvalN,
                             
                             permI = permI,
                             fixRankL = fixRankL,
                             
                             seedI = seedI,
                             plotSubC = plotSubC,
                             fig.pdfC = fig.pdfC,                   
                             info.txtC = info.txtC)
              
              y_name.c <- y
              
            } else if (is.list(y)) { # this mode is only used internally for MultiAssayExperiment
              
              stopifnot(length(y) == 1)
              stopifnot(!is.null(names(y)))
              
              bsg <- biosign(x = t(SummarizedExperiment::assay(x)),
                             y = y[[1]],
                             methodVc = methodVc,
                             bootI = bootI,
                             pvalN = pvalN,
                             
                             permI = permI,
                             fixRankL = fixRankL,
                             
                             seedI = seedI,
                             plotSubC = plotSubC,
                             fig.pdfC = fig.pdfC,                   
                             info.txtC = info.txtC)
              
              y_name.c <- names(y)
              
            } else
              stop("'y' must be the name of a column of the colData slot of the 'SummarizedExperiment' object")
        
            method.vc <- methodVc
            if (length(methodVc) == 1 && methodVc == "all")
              method.vc <- c("plsda", "forest", "svm")
            
            x@metadata[["biosign"]][[paste0(make.names(y_name.c),
                                            "_",
                                            paste(gsub("randomforest", "forest", method.vc),
                                                  collapse = "."))]] <- bsg
            
            if (!all(is.na(bsg@accuracyMN["S", ]))) {
              
              tierMC <- bsg@tierMC
              
              fdaDF <- SummarizedExperiment::rowData(x)
              
              for (colI in 1:ncol(tierMC)) {
                
                tierVc <- ropls:::.genVec(x, "feature", "character")
                
                tierVc[rownames(tierMC)] <- tierMC[, colI]
                
                fdaDF[, paste(make.names(y_name.c),
                              "biosign",
                              gsub("randomforest", "forest", colnames(tierMC)[colI]),
                              sep = "_")] <- tierVc
                
              }
              
              SummarizedExperiment::rowData(x) <- fdaDF
              
            }
            
            return(invisible(x))
            
          })

#### biosign (MultiAssayExperiment) ####

#' @rdname biosign
#' @export
setMethod("biosign", signature(x = "MultiAssayExperiment"),
          function(x,
                   y,
                   methodVc = c("all", "plsda", "randomforest", "svm")[1],
                   bootI = 50,
                   pvalN = 0.05,
                   
                   permI = 1,
                   fixRankL = FALSE,
                   
                   seedI = 123,
                   plotSubC = "",
                   fig.pdfC = c("none", "interactive", "myfile.pdf")[2],                   
                   info.txtC = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!is.character(y))
              stop("'y' must be a character when the 'biosign' method is applied to an 'MultiAssayExperiment' instance")
            
            
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
            
            if (!(y %in% colnames(MultiAssayExperiment::colData(x)))) {
              stop("'y' must be the name of a column for the colData of the MultiAssayExperiment instance")
            } else {
              colData.DF <- MultiAssayExperiment::colData(x)
              y.fcvcn <- colData.DF[, y]
              names(y.fcvcn) <- rownames(colData.DF)
              method.vc <- methodVc
              if (length(methodVc) == 1 && methodVc == "all")
                method.vc <- c("plsda", "forest", "svm")
              type.c <- paste0(make.names(y), "_", paste(gsub("randomforest", "forest", method.vc), collapse = "."))
            }
            
            for (setC in names(x)) {
              
              if (info.txtC != "none")
                cat("\n\nSelecting the features for the '", setC, "' dataset:\n", sep = "")
              
              plotL <- TRUE
              
              set.se <- x[[setC]]
              
              y_set.fcvcn <- y.fcvcn[colnames(set.se)]
              
              y_set.ls <- list(y_set.fcvcn)
              names(y_set.ls) <- y
              
              set.se <- biosign(set.se,
                                y = y_set.ls,
                                methodVc = methodVc,
                                bootI = bootI,
                                pvalN = pvalN,
                                
                                permI = permI,
                                fixRankL = fixRankL,
                                
                                seedI = seedI,
                                plotSubC = plotSubC,
                                fig.pdfC = figPdfC,                   
                                info.txtC = infTxtC)
              
              set.biosign <- getBiosign(set.se)[[type.c]]
              
              if (all(is.na(set.biosign@accuracyMN["S", ]))) {
                
                plotL <- FALSE
                
              }
              
              if (fig.pdfC != "none" && plotL) {
                
                plot(x = set.biosign,
                     y = y,
                     typeC = "tier",
                     
                     plotSubC = paste0("[", setC, "]"),                  
                     fig.pdfC = fig.pdfC,
                     info.txtC = info.txtC)
                
                plot(x = set.biosign,
                     y = y,
                     typeC = "boxplot",
                     
                     plotSubC = paste0("[", setC, "]"),                  
                     fig.pdfC = fig.pdfC,
                     info.txtC = info.txtC)
                
              }
              
              x[[setC]] <- set.se
              
            }
            
            if (!(fig.pdfC %in% c("none", "interactive")))
              grDevices::dev.off()
            
            if (!(info.txtC %in% c("none", "interactive")))
              sink()
            
            return(invisible(x))
            
          })


#### biosign (ExpressionSet) ####

#' @rdname biosign
#' @export
setMethod("biosign", signature(x = "ExpressionSet"),
          function(x,
                   y,
                   methodVc = c("all", "plsda", "randomforest", "svm")[1],
                   bootI = 50,
                   pvalN = 0.05,
                   
                   permI = 1,
                   fixRankL = FALSE,
                   
                   seedI = 123,
                   plotSubC = "",
                   fig.pdfC = c("none", "interactive", "myfile.pdf")[2],                   
                   info.txtC = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(y %in% colnames(Biobase::pData(x)))) {
              stop("'y' must be the name of a column of the phenoData slot of the 'ExpressionSet' object", call. = FALSE)
            } else {
              rspFcVc <- Biobase::pData(x)[, y]
              bsg <- biosign(x = t(Biobase::exprs(x)),
                             y = rspFcVc,
                             methodVc = methodVc,
                             bootI = bootI,
                             pvalN = pvalN,
                             
                             permI = permI,
                             fixRankL = fixRankL,
                             
                             seedI = seedI,
                             plotSubC = plotSubC,
                             fig.pdfC = fig.pdfC,                   
                             info.txtC = info.txtC)
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


#### biosign (MultiDataSet) ####

#' @rdname biosign
#' @export
setMethod("biosign", signature(x = "MultiDataSet"),
          function(x,
                   y,
                   methodVc = c("all", "plsda", "randomforest", "svm")[1],
                   bootI = 50,
                   pvalN = 0.05,
                   
                   permI = 1,
                   fixRankL = FALSE,
                   
                   seedI = 123,
                   plotSubC = "",
                   fig.pdfC = c("none", "interactive", "myfile.pdf")[2],                   
                   info.txtC = c("none", "interactive", "myfile.txt")[2]) {
            
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
              
              setBiosign <- tryCatch(biosigner::biosign(x = x[[setC]],
                                                        y = y,
                                                        methodVc = methodVc,
                                                        bootI = bootI,
                                                        pvalN = pvalN,
                                                        
                                                        permI = permI,
                                                        fixRankL = fixRankL,
                                                        
                                                        seedI = seedI,
                                                        plotSubC = plotSubC,
                                                        fig.pdfC = figPdfC,                   
                                                        info.txtC = infTxtC),
                                     error = function(e) NULL)
              
              if (all(is.na(setBiosign@accuracyMN["S", ]))) {
                
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
                
                plot(x = setBiosign,
                     y = y,
                     typeC = "tier",
                     
                     plotSubC = paste0("[", setC, "]"),                  
                     fig.pdfC = fig.pdfC,
                     info.txtC = info.txtC)
                
                plot(x = setBiosign,
                     y = y,
                     typeC = "boxplot",
                     
                     plotSubC = paste0("[", setC, "]"),                  
                     fig.pdfC = fig.pdfC,
                     info.txtC = info.txtC)
                
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


####     getAccuracyMN    ####

#' @rdname getAccuracyMN
#' @export
setMethod("getAccuracyMN", "biosign",
          function(object) {
            return(object@accuracyMN)
          })


#### getBiosign (SummarizedExperiment) ####

#' @rdname getBiosign
#' @export
setMethod("getBiosign", "SummarizedExperiment",
          function(object) {
            return(object@metadata[["biosign"]])
          })


#### getBiosign (MultiAssayExperiment) ####

#' @rdname getBiosign
#' @export
setMethod("getBiosign", "MultiAssayExperiment",
          function(object) {
            set.vc <- names(object)
            biosign.ls <- vector(mode = "list", length = length(set.vc))
            names(biosign.ls) <- set.vc
            for (set.c in set.vc) {
              biosign.ls[[set.c]] <- object[[set.c]]@metadata[["biosign"]]
            }
            return(biosign.ls)
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


####    getSignatureLs    ####

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
          function(object, newdata, tierMaxC = "S")
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






