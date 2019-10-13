####    plot  (biosignMultiDataSet)  ####

#' Plot method for biosign signatures
#'
#' This function plots signatures obtained by \code{biosign}.
#'
#' @aliases plot.biosignMultiDataSet plot,biosignMultiDataSet-method
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
#' # Plotting the selected signatures
#' plot(nciBiosign)
#' @rdname plot
#' @export
setMethod("plot", signature(x = "biosignMultiDataSet"),
          function(x,
                   y,
                   fig.pdfC = c("none", "interactive", "myfile.pdf")[2],
                   info.txtC = c("none", "interactive", "myfile.txt")[2],
                   ...) {
            
            if (!(info.txtC %in% c("none", "interactive"))) {
              sink(info.txtC, append = TRUE)
            }
            
            infTxtC <- info.txtC
            if (infTxtC != "none")
              infTxtC <- "interactive"
            
            if (!(fig.pdfC %in% c("none", "interactive"))) {
              pdf(fig.pdfC)
            }
            
            figPdfC <- fig.pdfC
            if (figPdfC != "none")
              figPdfC <- "interactive"
            
            biosignLs <- x@biosignLs
            
            for (setI in 1:length(biosignLs)) {
              plot(x@biosignLs[[setI]], subC = paste0("[", names(biosignLs)[setI], "]"),
                   typeC = "tier",
                   fig.pdfC = figPdfC,
                   info.txtC = infTxtC,
                   ...)
              plot(x@biosignLs[[setI]], subC = paste0("[", names(biosignLs)[setI], "]"),
                   typeC = "boxplot",
                   fig.pdfC = figPdfC,
                   info.txtC = infTxtC,
                   ...)
            }
            
            if (!(fig.pdfC %in% c("none", "interactive")))
              dev.off()
            
            if (!(info.txtC %in% c("none", "interactive")))
              sink()
            
          })

####    plot (biosign)   ####

#' Plot method for 'biosign' signature objects
#'
#' Displays classifier tiers or individual boxplots from selected features
#'
#' @aliases plot.biosign plot,biosign-method
#' @param x An S4 object of class \code{biosign}, created by the \code{biosign}
#' function.
#' @param y Currently not used.
#' @param tierMaxC Character: Maximum level of tiers to display: Either 'S' and
#' 'A', (for boxplot), or also 'B', 'C', 'D', and 'E' (for tiers) by decreasing
#' number of selections
#' @param typeC Character: Plot type; either 'tier' [default] displaying the
#' comparison of signatures up to the selected 'tierMaxC' or 'boxplot' showing
#' the individual boxplots of the features selected by all the classifiers
#' @param plotSubC Character: Graphic subtitle
#' @param fig.pdfC Character: File name with '.pdf' extension for the figure;
#' if 'interactive' (default), figures will be displayed interactively; if 'none',
#' no figure will be generated
#' @param info.txtC Character: File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @param file.pdfC Character: deprecated; use the 'fig.pdfC' argument instead
#' @param .sinkC Character: deprecated; use the 'info.txtC' argument instead
#' @param ... Currently not used.
#' @return A plot is created on the current graphics device.
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
#' plot(diaSign, typeC = "boxplot")
#'
#' detach(diaplasma)
#'
#' @rdname plot
#' @export
setMethod("plot", signature(x = "biosign"),
          function(x,
                   y,
                   tierMaxC = "S",
                   typeC = c("tier", "boxplot")[1],
                   
                   plotSubC = NA,                  
                   fig.pdfC = c("none", "interactive", "myfile.pdf")[2],
                   info.txtC = c("none", "interactive", "myfile.txt")[2],
                   
                   file.pdfC = NULL,
                   .sinkC = NULL,
                   ...) {
            
            if (!is.null(file.pdfC)) {
              warning("'file.pdfC' argument is deprecated; use 'fig.pdfC' instead.",
                      call. = FALSE)
              fig.pdfC <- file.pdfC
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
            
            if (fig.pdfC == "none")
              stop("'fig.pdfC' cannot be set to 'none' in the 'plot' method.",
                   call. = FALSE)
            
            if (!(info.txtC %in% c("none", "interactive")))
              sink(info.txtC, append = TRUE)
            
            tierFullVc <- c("S", LETTERS[1:5])
            
            if (length(tierMaxC) != 1 || any(!(tierMaxC %in% tierFullVc))) {
              stop("'tierMaxC' argument must be either '",
                   paste(tierFullVc, collapse = "', '"),
                   "' for the 'tier' plot.", call. = FALSE)
            } else if (typeC == "boxplot" && any(!(tierMaxC %in% c("S", "A")))) {
              stop("'tierMaxC' argument must be either 'S' or 'A' for the 'boxplot'.",
                   call. = FALSE)
            } else
              tierVc <- tierFullVc[1:which(tierFullVc == tierMaxC)]
            
            switch(typeC,
                   
                   tier = {
                     
                     if (sum(x@tierMC %in% tierVc) == 0) {
                       
                       stop("No signature up to tier '",
                            tierMaxC,
                            "' to be plotted.",
                            call. = FALSE)
                       
                     }
                     
                   },
                   boxplot = {
                     
                     switch(tierMaxC,
                            S = {
                              sgnAllVc = x@signatureLs[["complete"]]
                              xSubMN <- x@xSubMN
                            },
                            A = {
                              sgnAllVc = x@AS[["signatureLs"]][["complete"]]
                              xSubMN <- x@AS[["xSubMN"]]
                            })
                     
                     if (length(sgnAllVc) == 0) {
                       stop("No signature up to tier '",
                            tail(tierVc, 1), "' to be plotted.", call. = FALSE)
                     }
                     
                   })
            
            
            if (fig.pdfC != "interactive")
              pdf(fig.pdfC)
            
            opar <- par(no.readonly = TRUE)
            
            switch(typeC,
                   
                   tier = {
                     
                     if (length(setdiff(tierVc, c(x@tierMC)))) {
                       tierNotFoundVc <- tierFullVc[tierFullVc %in% setdiff(tierVc, c(x@tierMC))]
                       warning("tierMC does not contain the following values: '",
                               paste(tierNotFoundVc, collapse = "', '"), "'",
                               call. = FALSE)
                     }
                     inputMC <- x@tierMC[apply(x@tierMC, 1,
                                               function(rowVc) sum(rowVc %in% tierVc) > 0), ,
                                         drop = FALSE]
                     tierFullVi <- 1:length(tierFullVc)
                     names(tierFullVi) <- tierFullVc
                     inputMN <- inputMC
                     for (j in 1:ncol(inputMN))
                       inputMN[, j] <- tierFullVi[inputMC[, j]]
                     mode(inputMN) <- "numeric"
                     
                     imageMN <- inputMN
                     
                     colnames(imageMN) <- gsub("plsda", "PLS-DA",
                                               gsub("randomforest", "Random Forest",
                                                    gsub("svm", "SVM", colnames(imageMN))))
                     
                     if (length(rownames(imageMN)) == 0) {
                       
                       rownameInVc <- rep("", times = nrow(imageMN))
                       
                     } else {
                       
                       rownameInVc <- rownames(imageMN)
                       
                       warnOpN <- getOption("warn")
                       options(warn = -1)
                       rownameInCharVsNumL <- NA %in% as.numeric(rownameInVc)
                       options(warn = warnOpN)
                       
                       if (rownameInCharVsNumL) {
                         
                         rownameInDuplicateVc <- duplicated(rownameInVc)
                         
                         rownameInVc <- sapply(1:length(rownameInVc), function(i) ifelse(!rownameInDuplicateVc[i],
                                                                                         yes = rownameInVc[i],
                                                                                         no = ""))
                         
                       }
                     }
                     
                     colnameInVc <- colnames(imageMN)
                     
                     
                     imageMN <- t(imageMN[nrow(imageMN):1, , drop = FALSE])
                     
                     par(bg = "white",
                         font = 2,
                         lwd  = 2)
                     
                     layout(matrix(c(2, 1),
                                   nrow = 1,
                                   ncol = 2,
                                   byrow = TRUE),
                            widths = c(8, 2))
                     
                     tierFullColVc <- c("#1A9850", "#91CF60", "#D9EF8B", "#FEE08B", "#FC8D59", "#D73027")
                     names(tierFullColVc) <- tierFullVc
                     
                     ## draw color scale
                     
                     par(mar = c(1.1, 0.6, 7.6, 4.1))
                     
                     plot(x = 0,
                          y = 0,
                          font.axis = 2,
                          font.lab = 2,
                          type = "n",
                          xlim = c(0, 1),
                          ylim = c(0, 6),
                          xlab = "",
                          ylab = "",
                          xaxs = "i",
                          yaxs = "i",
                          xaxt = "n",
                          yaxt = "n")
                     
                     rect(xleft = 0,
                          ybottom = 0:5,
                          xright = 1,
                          ytop = 1:6,
                          col = rev(tierFullColVc),
                          border = NA)
                     
                     axis(at = tierFullVi - 0.5,
                          font = 2,
                          font.axis = 2,
                          labels = rev(tierFullVc),
                          las = 1,
                          lwd = 2,
                          lwd.ticks = 2,
                          side = 4)
                     
                     arrows(par("usr")[2],
                            par("usr")[4],
                            par("usr")[2],
                            par("usr")[3],
                            code = 0,
                            lwd = 2,
                            xpd = TRUE)
                     
                     ## draw image
                     
                     par(mar = c(1.1,
                                 17,
                                 7.6,
                                 0.3))
                     
                     image(x = 1:nrow(imageMN),
                           y = 1:ncol(imageMN),
                           z = imageMN,
                           col = tierFullColVc[min(inputMN):max(inputMN)],
                           font.axis = 2,
                           font.lab = 2,
                           xaxt = "n",
                           yaxt = "n",
                           xlab = "",
                           ylab = "")
                     
                     xLabelAtVn <- 1:ncol(imageMN)
                     yLabelAtVn <- 1:nrow(imageMN)
                     
                     ## if(colnameInCharVsNumL) {
                     
                     xLabelVc <- colnameInVc[colnameInVc != ""]
                     
                     xLabelLowIndiceVn <- which(colnameInVc != "")
                     
                     xLabelSpanVn <- diff(c(xLabelLowIndiceVn, length(colnameInVc) + 1))
                     
                     xLabelAtVn <- xLabelLowIndiceVn - rep(1, times = length(xLabelLowIndiceVn)) + xLabelSpanVn / 2 + rep(0.5, times = length(xLabelLowIndiceVn))
                     
                     par(cex = 1)
                     
                     axis(side = 3,
                          at = xLabelAtVn,
                          font = 2,
                          labels = xLabelVc,
                          las = 2,
                          line = -0.5,
                          tick = FALSE)
                     
                     par(cex = 1)
                     
                     
                     if (rownameInCharVsNumL) {
                       
                       yLabelVc <- rownameInVc[rownameInVc != ""]
                       
                       yLabelLowIndiceVn <- which(rownameInVc != "")
                       
                       yLabelSpanVn <- diff(c(yLabelLowIndiceVn, length(rownameInVc) + 1))
                       
                       yLabelAtVn <- ncol(imageMN) - rev(yLabelLowIndiceVn) + rep(1, times = length(yLabelLowIndiceVn)) - rev(yLabelSpanVn) / 2 + rep(0.5, times = length(yLabelLowIndiceVn))
                       
                       par(cex = 1)
                       
                       axis(side = 2,
                            at = yLabelAtVn,
                            font = 2,
                            hadj = 1,
                            labels = rev(yLabelVc),
                            las = 2,
                            line = -0.5,
                            tick = FALSE)
                       
                     } else {
                       
                       rownameVn <- as.numeric(rownameInVc)
                       
                       prettyVn <- pretty(rownameVn)
                       
                       prettyVn <- prettyVn[min(rownameVn) <= prettyVn & prettyVn <= max(rownameVn)]
                       
                       indiceVn <- numeric()
                       
                       for (k in 1:length(prettyVn))
                         indiceVn[k] <- which(abs(rownameVn - prettyVn[k]) == min(abs(rownameVn - prettyVn[k])))[1]
                       
                       axis(side = 2,
                            at = nrow(imageMN) - rev(indiceVn),
                            font = 2,
                            labels = as.character(prettyVn))
                       
                     }
                     
                     par(cex = 1)
                     
                     ## additional lines
                     
                     if (rownameInCharVsNumL)
                       abline(h = ncol(imageMN) - which(rownameInVc != "") + 1 + 0.5,
                              lwd = 1)
                     
                     ## if (colnameInCharVsNumL)
                     abline(v = which(colnameInVc != "") - 0.5,
                            lwd = 1)
                     
                     
                     ## border
                     
                     box(lwd = 2)
                     
                     ## arrows at the end of the axes
                     
                     arrows(par("usr")[1],
                            par("usr")[4],
                            par("usr")[1],
                            par("usr")[3],
                            length = 0.1,
                            lwd = 2,
                            xpd = TRUE)
                     
                     box(lwd = 3)
                     
                     ## Title
                     
                     mainC <- "Tiers of the selected features"
                     if (!is.na(plotSubC))
                       mainC <- paste0(mainC, " \n", plotSubC)
                     
                     title(mainC, adj = 0.1, outer = TRUE, line = -2.5)
                     
                   },
                   boxplot = {
                     
                     palVc <- c(code100_plsda = "#CB181D",
                                ## code110_plsda_randomforest = "#F16913",
                                code010_randomforest = "#238B45",
                                ## code011_randomforest_svm = "#807DBA",
                                code001_svm = "#2171B5")
                     ## code101_plsda_svm = "          #BF812D",
                     ## code111_plsda_randomforest_svm = "black"
                     ## )
                     
                     palDF <- data.frame(row.names = substr(names(palVc), 5, 7),
                                         classifier = substr(names(palVc), 9, nchar(names(palVc))),
                                         color = palVc,
                                         stringsAsFactors = FALSE)
                     
                     dimN <- ceiling(sqrt(length(sgnAllVc)))
                     layout(matrix(1:dimN^2, nrow = dimN, byrow = TRUE))
                     par(font = 2, font.axis = 2, font.lab = 2, las = 1,
                         mar = c(2.1, 2.6, 2.6, 1.1),
                         oma = c(0, 0, 2.1, 0))
                     for (bmkC in sgnAllVc) {
                       bmkTirVc <- x@tierMC[bmkC, ]
                       names(bmkTirVc) <- colnames(x@tierMC)
                       bmkTirSigVl <- bmkTirVc %in% tierVc
                       bmkModVc <- paste(names(bmkTirVc)[bmkTirSigVl], collapse = ", ")
                       bmkModVc <- gsub("plsda", "PLSDA", gsub("randomforest", 'RF', gsub("svm", "SVM", bmkModVc)))
                       bmkColC <- ifelse(sum(bmkTirSigVl) == 1,
                                         palDF[paste(as.numeric(bmkTirSigVl), collapse = ""), "color"],
                                         "black")
                       boxplot(xSubMN[, bmkC] ~ x@yFc,
                               border = bmkColC,
                               main = "")
                       mtext(ifelse(nchar(bmkC) > 23, paste0(substr(bmkC, 1, 20), "."), bmkC),
                             cex = 0.8,
                             line = 1.1)
                       mtext(bmkModVc, line = 0.1, cex = 0.6)
                     }
                     mainC <- paste0("'", switch(tierMaxC, S = "S", A = "S+A"), "' signature")
                     if (!is.na(plotSubC))
                       mainC <- paste0(mainC, " ", plotSubC)
                     title(main = mainC,
                           line = 0.5, cex.main = 1.5, outer = TRUE)
                     
                   })
            
            par(opar)
            
            if (fig.pdfC != "interactive")
              dev.off()
            
            
            ## Closing connection
            
            if (!(info.txtC %in% c("none", "interactive")))
              sink()
            
          })
