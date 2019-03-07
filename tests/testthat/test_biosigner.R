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
                       fig.pdfC = NULL)
  
  set.seed(NULL)
  
  plot(biosignLs, fig.pdfC = "test.pdf")
  plot(biosignLs, tierMaxC = "A", fig.pdfC = "test.pdf")
  plot(biosignLs, typeC = "boxplot", fig.pdfC = "test.pdf")
  plot(biosignLs, tierMaxC = "A", typeC = "boxplot", fig.pdfC = "test.pdf")
  
  if (.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    testthat::expect_identical(biosignLs@tierMC["m427.215t07.9", "plsda"],
                               "A")
    testthat::expect_equivalent(biosignLs@accuracyMN["S", "plsda"],
                                0.7365702, tolerance = 1e-7)
    
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
                                0.7393993, tolerance = 1e-7)
    testthat::expect_equivalent(biosignLs@modelLs[["randomforest"]][["votes"]][2],
                                0.005988024, tolerance = 1e-6)
    
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
                               "A")
    testthat::expect_equivalent(biosignLs@accuracyMN["AS", "svm"],
                                0.8313882, tolerance = 1e-7)
    testthat::expect_equivalent(biosignLs@AS[["modelLs"]][["svm"]][["rho"]],
                                -0.7013267, tolerance = 1e-6)
    
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
    
    accMC <- matrix(c("0.7234427", "0.6964646", "0.7620893", "0.7015924", "0.7970713", "0.714572", "0.7458978", "0.5143775", NA),
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
  
  
  sigLs <- list(plsda = c("Oxoglutaric acid", "p-Anisic acid", "Testosterone glucuronide", "Acetylphenylalanine", "Malic acid", "Pantothenic acid", "Gluconic acid and/or isomers", "alpha-N-Phenylacetyl-glutamine", "Citric acid", "Glucuronic acid and/or isomers", "Hippuric acid", "Phe-Tyr-Asp (and isomers)", "Threonic acid/Erythronic acid"),
                randomforest = c("Oxoglutaric acid", "p-Anisic acid", "Testosterone glucuronide"),
                svm = c("Oxoglutaric acid", "p-Anisic acid", "Testosterone glucuronide", "Acetylphenylalanine", "Malic acid", "Pantothenic acid", "Taurine", "N4-Acetylcytidine", "Monoethyl phthalate"),
                complete = c("Oxoglutaric acid", "p-Anisic acid", "Testosterone glucuronide", "Acetylphenylalanine", "Malic acid", "Pantothenic acid", "Gluconic acid and/or isomers", "Taurine", "N4-Acetylcytidine", "alpha-N-Phenylacetyl-glutamine", "Citric acid", "Glucuronic acid and/or isomers", "Hippuric acid", "Monoethyl phthalate", "Phe-Tyr-Asp (and isomers)", "Threonic acid/Erythronic acid"))
  
  if (.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    testthat::expect_identical(getSignatureLs(biosignLs), sigLs)
    
    accMC <- matrix(c("0.8958831", "0.8622427", "0.8690036",
                      "0.8933625", "0.8792595", "0.8688045",
                      "0.8974811", "0.9100602", "0.8964705"),
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


