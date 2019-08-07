library(biosigner)

context("Testing 'biosigner'")

test_that("biosign_plsda", {
  
  data(diaplasma)
  
  varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))
  
  set.seed(123)
  
  biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                       diaplasma[["sampleMetadata"]][, "type"],
                       methodVc = "plsda",
                       bootI = 5,
                       fig.pdfC = "none")
  
  set.seed(NULL)
  
  plot(biosignLs, fig.pdfC = "test.pdf")
  plot(biosignLs, tierMaxC = "A", fig.pdfC = "test.pdf")
  plot(biosignLs, typeC = "boxplot", fig.pdfC = "test.pdf")
  plot(biosignLs, tierMaxC = "A", typeC = "boxplot", fig.pdfC = "test.pdf")
  
  if (.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    testthat::expect_identical(biosignLs@tierMC["m427.215t07.9", "plsda"],
                               "A")
    testthat::expect_equivalent(biosignLs@accuracyMN["S", "plsda"],
                                0.6520116, tolerance = 1e-7)
    
    library(ropls)
    
    testthat::expect_equivalent(getSummaryDF(biosignLs@modelLs[["plsda"]])[, "Q2(cum)"],
                                0.271, tolerance = 1e-6)
    
  }
  
})

test_that("biosign_randomforest", {
  
  data(diaplasma)
  
  varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))
  
  set.seed(123)
  
  biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                       diaplasma[["sampleMetadata"]][, "type"],
                       methodVc = "randomforest",
                       bootI = 5,
                       fig.pdfC = "test.pdf")
  
  set.seed(NULL)
  
  if (.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    
    testthat::expect_identical(biosignLs@tierMC["m427.215t07.9", "randomforest"],
                               "S")
    
    testthat::expect_equivalent(biosignLs@accuracyMN["AS", "randomforest"],
                                0.734806, tolerance = 1e-7)
    testthat::expect_equivalent(biosignLs@modelLs[["randomforest"]][["votes"]][2],
                                0, tolerance = 1e-6)
    
  }
  
})

test_that("biosign_svm", {
  
  data(diaplasma)
  
  varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))
  
  set.seed(123)
  
  biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                       diaplasma[["sampleMetadata"]][, "type"],
                       methodVc = "svm",
                       bootI = 5,
                       fig.pdfC = "test.pdf")
  
  set.seed(NULL)
  
  if (.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    testthat::expect_identical(biosignLs@tierMC["m123.998t01.0", "svm"],
                               "E")
    testthat::expect_equivalent(biosignLs@accuracyMN["AS", "svm"],
                                0.7443567, tolerance = 1e-7)
    testthat::expect_equivalent(biosignLs@AS[["modelLs"]][["svm"]][["rho"]],
                                -0.7788987, tolerance = 1e-6)
    
  }
})

test_that("biosign_predict", {
  
  data(diaplasma)
  
  varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))
  
  samTotI <- nrow(diaplasma[["dataMatrix"]])
  trainVi <- 1:floor(samTotI/2)
  
  set.seed(123)
  
  biosignLs <- biosign(diaplasma[["dataMatrix"]][trainVi, varSelVi],
                       diaplasma[["sampleMetadata"]][trainVi, "type"],
                       bootI = 1,
                       fig.pdfC = "test.pdf")
  
  set.seed(NULL)
  
  predDF <- predict(biosignLs,
                    diaplasma[["dataMatrix"]][setdiff(1:samTotI, trainVi), varSelVi])
  predDIA043Vc <- as.character(unlist(predDF["DIA043", ]))
  
  if (.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    testthat::expect_identical(predDIA043Vc,
                               c("T2", "T1", "T2"))
    
  }
})

test_that("biosign_diaplasma", {
  
  data(diaplasma)
  
  varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))
  
  set.seed(123)
  
  biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                       diaplasma[["sampleMetadata"]][, "type"],
                       bootI = 5,
                       fig.pdfC = "test.pdf")
  
  set.seed(NULL)
  
  sigLs <- list(plsda = "m189.040t01.2",
                randomforest = "m427.215t07.9",
                svm = character(0),
                complete = c("m427.215t07.9", "m189.040t01.2"))
  
  if (.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    testthat::expect_identical(getSignatureLs(biosignLs), sigLs)
    
    accMC <- matrix(c("0.7178939", "0.7464452", "0.737702", "0.703651", "0.8389099", "0.6439541", "0.7066359", "0.6589018", NA),
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
  
  library(ropls)
  
  data(sacurine)
  
  set.seed(123)
  
  biosignLs <- biosign(sacurine[["dataMatrix"]],
                       sacurine[["sampleMetadata"]][, "gender"],
                       bootI = 5,
                       fig.pdfC = "test.pdf")
  
  set.seed(NULL)
  
  
  sigLs <- list(plsda = c("p-Anisic acid", "Testosterone glucuronide", "Malic acid", "Pantothenic acid"),
                randomforest = c("p-Anisic acid", "Testosterone glucuronide", "Oxoglutaric acid"),
                svm = c("p-Anisic acid", "Testosterone glucuronide", "Oxoglutaric acid", "N2-Acetylaminoadipic acid", "N4-Acetylcytidine"),
                complete = c("p-Anisic acid", "Testosterone glucuronide", "Oxoglutaric acid", "Malic acid", "N2-Acetylaminoadipic acid", "N4-Acetylcytidine", "Pantothenic acid"))
  
  if (.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    testthat::expect_identical(getSignatureLs(biosignLs), sigLs)
    
    accMC <- matrix(c("0.8713552", "0.9211478", "0.8373789", "0.8539239", "0.8675991", "0.8644966", "0.9086947", "0.8895256", "0.8724575"),
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
  attach(diaplasma)

  diaSet <- Biobase::ExpressionSet(assayData = t(dataMatrix),
                                   phenoData = new("AnnotatedDataFrame",
                                                   data = sampleMetadata),
                                   featureData = new("AnnotatedDataFrame",
                                                     data = variableMetadata),
                                   experimentData = new("MIAME",
                                                        title = "diaplasma"))
  
  featureSelVl <- variableMetadata[, "mzmed"] >= 490 & variableMetadata[, "mzmed"] < 500
  diaSet <- diaSet[featureSelVl, ]
  
  set.seed(123)
  diaSign <- biosign(diaSet, "type", bootI = 5, fig.pdfC = "test.pdf")
  
  diaSet <- biosigner::getEset(diaSign)
  
  testthat::expect_equivalent(Biobase::fData(diaSet)["m490.111t12.5", "rtmed"],
                              749.17605, tolerance = 1e-5)
  
  testthat::expect_identical(Biobase::fData(diaSet)["m490.111t12.5", "type_biosign_forest"],
                             "E")
  
})


