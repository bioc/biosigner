library(biosigner)

context("Testing 'biosigner'")

test_that("biosign_plsda", {
  
  data(diaplasma)
  
  varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))
  
  biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                       diaplasma[["sampleMetadata"]][, "type"],
                       methodVc = "plsda",
                       bootI = 5,
                       fig.pdfC = "none")
  
  plot(biosignLs, fig.pdfC = "test.pdf")
  plot(biosignLs, tierMaxC = "A", fig.pdfC = "test.pdf")
  plot(biosignLs, typeC = "boxplot", fig.pdfC = "test.pdf")
  plot(biosignLs, tierMaxC = "A", typeC = "boxplot", fig.pdfC = "test.pdf")
  
  if (.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    testthat::expect_identical(biosignLs@tierMC["m427.215t07.9", "plsda"],
                               "A")
    testthat::expect_equivalent(biosignLs@accuracyMN["S", "plsda"],
                                0.6968404, tolerance = 1e-7)
    testthat::expect_equivalent(getSummaryDF(biosignLs@modelLs[["plsda"]])[, "Q2(cum)"],
                                0.271, tolerance = 1e-6)
    
  }
  
})

test_that("biosign_randomforest", {
  
  data(diaplasma)
  
  varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))
  
  biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                       diaplasma[["sampleMetadata"]][, "type"],
                       methodVc = "randomforest",
                       bootI = 5,
                       fig.pdfC = "test.pdf")
  
  if (.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    testthat::expect_identical(biosignLs@tierMC["m427.215t07.9", "randomforest"],
                               "S")
    testthat::expect_equivalent(biosignLs@accuracyMN["AS", "randomforest"],
                                0.7140211, tolerance = 1e-7)
    testthat::expect_equivalent(biosignLs@modelLs[["randomforest"]][["votes"]][2],
                                0, tolerance = 1e-6)
    
  }
  
})

test_that("biosign_svm", {
  
  data(diaplasma)
  
  varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))
  
  biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                       diaplasma[["sampleMetadata"]][, "type"],
                       methodVc = "svm",
                       bootI = 5,
                       fig.pdfC = "test.pdf")
  
  if (.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    testthat::expect_identical(biosignLs@tierMC["m123.998t01.0", "svm"],
                               "E")
    testthat::expect_equivalent(biosignLs@accuracyMN["AS", "svm"],
                                0.7447254, tolerance = 1e-7)
    testthat::expect_equivalent(biosignLs@AS[["modelLs"]][["svm"]][["rho"]],
                                -0.7788987, tolerance = 1e-6)
    
  }
})

test_that("biosign_predict", {
  
  data(diaplasma)
  
  varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))
  
  samTotI <- nrow(diaplasma[["dataMatrix"]])
  trainVi <- 1:floor(samTotI/2)
  
  biosignLs <- biosign(diaplasma[["dataMatrix"]][trainVi, varSelVi],
                       diaplasma[["sampleMetadata"]][trainVi, "type"],
                       bootI = 1,
                       fig.pdfC = "test.pdf")
  
  predDF <- predict(biosignLs,
                    diaplasma[["dataMatrix"]][setdiff(1:samTotI, trainVi), varSelVi])
  predDIA043Vc <- as.character(unlist(predDF["DIA043", ]))
  
  if (.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    testthat::expect_identical(predDIA043Vc,
                               c("T2", "T2", "T2"))
    
  }
})

test_that("biosign_diaplasma", {
  
  data(diaplasma)
  
  varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))
  
  biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                       diaplasma[["sampleMetadata"]][, "type"],
                       bootI = 5,
                       fig.pdfC = "test.pdf")
  
  sigLs <- list(plsda = "m189.040t01.2",
                randomforest = "m427.215t07.9",
                svm = character(0),
                complete = c("m427.215t07.9", "m189.040t01.2"))
  
  if (.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    testthat::expect_identical(getSignatureLs(biosignLs), sigLs)
    
    accMC <- matrix(c("0.7178939", "0.711159", "0.6968404", "0.7432468", "0.7140211", "0.6650198", "0.6881863", "0.7447254", NA),
                    nrow = 3,
                    ncol = 3,
                    dimnames = list(c("Full", "AS", "S"),
                                    c("plsda", "randomforest", "svm")))
    
    biosignMC <- round(getAccuracyMN(biosignLs), 7)
    mode(biosignMC) <- "character"
    
    testthat::expect_identical(biosignMC, accMC)
    
  }
  
})

test_that("biosign_sacurine", {
  
  data(sacurine, package = "ropls")
  
  biosignLs <- biosign(sacurine[["dataMatrix"]],
                       sacurine[["sampleMetadata"]][, "gender"],
                       bootI = 5,
                       fig.pdfC = "test.pdf")
  
  sigLs <- list(plsda = c("p-Anisic acid", "Testosterone glucuronide", "Pantothenic acid", "Malic acid"),
                randomforest = c("p-Anisic acid", "Testosterone glucuronide", "Oxoglutaric acid"),
                svm = c("p-Anisic acid", "Testosterone glucuronide", "Pantothenic acid", "Oxoglutaric acid", "Glucuronic acid and/or isomers", "2-Isopropylmalic acid", "3,7-Dimethyluric acid", "4-Acetamidobutanoic acid isomer 3", "N-Acetylleucine", "N2-Acetylaminoadipic acid", "N4-Acetylcytidine", "Pyrroledicarboxylic acid", "Taurine", "Xanthosine"),
                complete = c("p-Anisic acid", "Testosterone glucuronide", "Pantothenic acid", "Malic acid", "Oxoglutaric acid", "Glucuronic acid and/or isomers", "2-Isopropylmalic acid", "3,7-Dimethyluric acid", "4-Acetamidobutanoic acid isomer 3", "N-Acetylleucine", "N2-Acetylaminoadipic acid", "N4-Acetylcytidine", "Pyrroledicarboxylic acid", "Taurine", "Xanthosine"))
  
  if (.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    testthat::expect_identical(getSignatureLs(biosignLs), sigLs)
    
    accMC <- matrix(c("0.8713552", "0.8825936", "0.8768834", "0.8529397", "0.9200604", "0.8551306", "0.8842949", "0.927461", "0.9306658"),
                    nrow = 3,
                    ncol = 3,
                    dimnames = list(c("Full", "AS", "S"),
                                    c("plsda", "randomforest", "svm")))
    
    biosignMC <- round(getAccuracyMN(biosignLs), 7)
    mode(biosignMC) <- "character"
    
    testthat::expect_identical(biosignMC, accMC)
    
  }
  
})

test_that("ExpressionSet", {
  
  data(diaplasma)

  diaSet <- Biobase::ExpressionSet(assayData = t(diaplasma[["dataMatrix"]]),
                                   phenoData = new("AnnotatedDataFrame",
                                                   data = diaplasma[["sampleMetadata"]]),
                                   featureData = new("AnnotatedDataFrame",
                                                     data = diaplasma[["variableMetadata"]]),
                                   experimentData = new("MIAME",
                                                        title = "diaplasma"))
  
  featureSelVl <- diaplasma[["variableMetadata"]][, "mzmed"] >= 490 & diaplasma[["variableMetadata"]][, "mzmed"] < 500
  diaSet <- diaSet[featureSelVl, ]
  
  diaSign <- biosign(diaSet, "type", bootI = 5, fig.pdfC = "test.pdf")
  
  diaSet <- biosigner::getEset(diaSign)
  
  testthat::expect_equivalent(Biobase::fData(diaSet)["m490.111t12.5", "rtmed"],
                              749.17605, tolerance = 1e-5)
  
  testthat::expect_identical(Biobase::fData(diaSet)["m490.111t12.5", "type_biosign_forest"],
                             "E")
  
})

test_that("MultiDataSet", {
  
  # Loading the 'NCI60_4arrays' from the 'omicade4' package
  data("NCI60_4arrays", package = "omicade4")
  # Selecting two of the four datasets
  setNamesVc <- c("agilent", "hgu95")
  # Creating the MultiDataSet instance
  nciMset <- MultiDataSet::createMultiDataSet()
  # Adding the two datasets as ExpressionSet instances
  for (setC in setNamesVc) {
    # Getting the data
    exprMN <- as.matrix(NCI60_4arrays[[setC]])
    pdataDF <- data.frame(row.names = colnames(exprMN),
                          cancer = substr(colnames(exprMN), 1, 2),
                          stringsAsFactors = FALSE)
    fdataDF <- data.frame(row.names = rownames(exprMN),
                          name = rownames(exprMN),
                          stringsAsFactors = FALSE)
    # Building the ExpressionSet
    eset <- Biobase::ExpressionSet(assayData = exprMN,
                                   phenoData = new("AnnotatedDataFrame",
                                                   data = pdataDF),
                                   featureData = new("AnnotatedDataFrame",
                                                     data = fdataDF),
                                   experimentData = new("MIAME",
                                                        title = setC))
    # Adding to the MultiDataSet
    nciMset <- MultiDataSet::add_eset(nciMset, eset, dataset.type = setC,
                                      GRanges = NA, warnings = FALSE)
  }
  # Restricting to the 'ME' and 'LE' cancer types
  sampleNamesVc <- Biobase::sampleNames(nciMset[["agilent"]])
  cancerTypeVc <- Biobase::pData(nciMset[["agilent"]])[, "cancer"]
  nciMset <- nciMset[sampleNamesVc[cancerTypeVc %in% c("ME", "LE")], ]
  # Summary of the MultiDataSet
  nciMset
  # Selecting the significant features for PLS-DA, RF, and SVM classifiers, and getting back the updated MultiDataSet
  nciBiosign <- biosigner::biosign(nciMset, "cancer")
  nciMset <- biosigner::getMset(nciBiosign)
  # In the updated MultiDataSet, the updated featureData now contains the cancer_biosign_'classifier' columns
  # indicating the selected features
  # lapply(Biobase::fData(nciMset), head)
  
  testthat::expect_identical(Biobase::fData(nciMset)[["agilent"]]["GTPBP5", "cancer_biosign_plsda"],
                             "E")
  testthat::expect_identical(Biobase::fData(nciMset)[["hgu95"]]["TSPAN4", "cancer_biosign_plsda"],
                             "A")
  
  
})


