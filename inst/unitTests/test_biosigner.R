test_biosign_plsda <- function(){
  
  data(diaplasma)
  
  varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))
  
  set.seed(123)
  
  biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                       diaplasma[["sampleMetadata"]][, "type"],
                       methodVc = "plsda",
                       bootI = 5,
                       plotL = FALSE)
  
  set.seed(NULL)
  
  # plot(biosignLs)
  # plot(biosignLs, tierMaxC = "A")
  # plot(biosignLs, typeC = "boxplot")
  # plot(biosignLs, tierMaxC = "A", typeC = "boxplot")
  
  if(.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    checkEquals(biosignLs@tierMC["m427.215t07.9", "plsda"],
                "A")
    checkEqualsNumeric(biosignLs@accuracyMN["S", "plsda"],
                       0.7365702, tolerance = 1e-7)
    
    library(ropls)
    
    checkEqualsNumeric(getSummaryDF(biosignLs@modelLs[["plsda"]])[, "Q2(cum)"],
                       0.271, tolerance = 1e-6)
    
  }
  
}

test_biosign_randomforest <- function(){
  
  data(diaplasma)
  
  varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))
  
  set.seed(123)
  
  biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                       diaplasma[["sampleMetadata"]][, "type"],
                       methodVc = "randomforest",
                       bootI = 5,
                       plotL = FALSE)
  
  set.seed(NULL)
  
  if(.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    
    checkEquals(biosignLs@tierMC["m427.215t07.9", "randomforest"],
                "S")
    
    checkEqualsNumeric(biosignLs@accuracyMN["AS", "randomforest"],
                       0.7393993, tolerance = 1e-7)
    checkEqualsNumeric(biosignLs@modelLs[["randomforest"]][["votes"]][2],
                       0.005988024, tolerance = 1e-6)
    
  }
  
}

test_biosign_svm <- function(){
  
  data(diaplasma)
  
  varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))
  
  set.seed(123)
  
  biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                       diaplasma[["sampleMetadata"]][, "type"],
                       methodVc = "svm",
                       bootI = 5,
                       plotL = FALSE)
  
  set.seed(NULL)
  
  if(.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    checkEquals(biosignLs@tierMC["m123.998t01.0", "svm"],
                "A")
    checkEqualsNumeric(biosignLs@accuracyMN["AS", "svm"],
                       0.8313882, tolerance = 1e-7)
    checkEqualsNumeric(biosignLs@AS[["modelLs"]][["svm"]][["rho"]],
                       -0.7013267, tolerance = 1e-6)
    
  }
}

test_biosign_predict <- function() {
  
  data(diaplasma)
  
  varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))
  
  samTotI <- nrow(diaplasma[["dataMatrix"]])
  trainVi <- 1:floor(samTotI/2)
  
  set.seed(123)
  
  biosignLs <- biosign(diaplasma[["dataMatrix"]][trainVi, varSelVi],
                       diaplasma[["sampleMetadata"]][trainVi, "type"],
                       bootI = 1,
                       plotL = FALSE)
  
  set.seed(NULL)
  
  predDF <- predict(biosignLs,
                    diaplasma[["dataMatrix"]][setdiff(1:samTotI, trainVi), varSelVi])
  predDIA043Vc <- as.character(unlist(predDF["DIA043", ]))
  
  if(.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    checkEquals(predDIA043Vc,
                c("T2", "T1", "T2"))
    
  }
}

test_biosign_diaplasma <- function() {
  
  data(diaplasma)
  
  varSelVi <- seq(1, ncol(diaplasma[["dataMatrix"]]), by = round(ncol(diaplasma[["dataMatrix"]]) / 100))
  
  set.seed(123)
  
  biosignLs <- biosign(diaplasma[["dataMatrix"]][, varSelVi],
                       diaplasma[["sampleMetadata"]][, "type"],
                       bootI = 5,
                       plotL = FALSE)
  
  set.seed(NULL)
  
  sigLs <- list(plsda = "m189.040t01.2",
                randomforest = "m427.215t07.9",
                svm = character(0),
                complete = c("m427.215t07.9", "m189.040t01.2"))
  
  if(.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    checkIdentical(getSignatureLs(biosignLs), sigLs)
    
    accMC <- matrix(c("0.7234427", "0.6964646", "0.7620893", "0.7015924", "0.7970713", "0.714572", "0.7458978", "0.5143775", NA),
                    nrow = 3,
                    ncol = 3,
                    dimnames = list(c("Full", "AS", "S"),
                                    c("plsda", "randomforest", "svm")))
    
    biosignMC <- round(getAccuracyMN(biosignLs), 7)
    mode(biosignMC) <- "character"
    
    checkIdentical(biosignMC, accMC)
    
  }
  
}

test_biosign_sacurine <- function() {
  
  library(ropls)
  
  data(sacurine)
  
  set.seed(123)
  
  biosignLs <- biosign(sacurine[["dataMatrix"]],
                       sacurine[["sampleMetadata"]][, "gender"],
                       bootI = 5,
                       plotL = FALSE)
  
  set.seed(NULL)
  
  
  sigLs <- list(plsda = c("Oxoglutaric acid", "p-Anisic acid", "Testosterone glucuronide", "Acetylphenylalanine", "Malic acid", "Pantothenic acid", "Gluconic acid and/or isomers", "alpha-N-Phenylacetyl-glutamine", "Citric acid", "Glucuronic acid and/or isomers", "Hippuric acid", "Phe-Tyr-Asp (and isomers)", "Threonic acid/Erythronic acid"),
                randomforest = c("Oxoglutaric acid", "p-Anisic acid", "Testosterone glucuronide"),
                svm = c("Oxoglutaric acid", "p-Anisic acid", "Testosterone glucuronide", "Acetylphenylalanine", "Malic acid", "Pantothenic acid", "Taurine", "N4-Acetylcytidine", "Monoethyl phthalate"),
                complete = c("Oxoglutaric acid", "p-Anisic acid", "Testosterone glucuronide", "Acetylphenylalanine", "Malic acid", "Pantothenic acid", "Gluconic acid and/or isomers", "Taurine", "N4-Acetylcytidine", "alpha-N-Phenylacetyl-glutamine", "Citric acid", "Glucuronic acid and/or isomers", "Hippuric acid", "Monoethyl phthalate", "Phe-Tyr-Asp (and isomers)", "Threonic acid/Erythronic acid"))
  
  if(.Platform$OS.type != "windows" || .Machine$sizeof.pointer == 8) {
    
    checkIdentical(getSignatureLs(biosignLs), sigLs)
    
    accMC <- matrix(c("0.8958831", "0.8622427", "0.8690036",
                      "0.8933625", "0.8792595", "0.8688045",
                      "0.8974811", "0.9100602", "0.8964705"),
                    nrow = 3,
                    ncol = 3,
                    dimnames = list(c("Full", "AS", "S"),
                                    c("plsda", "randomforest", "svm")))
    
    biosignMC <- round(getAccuracyMN(biosignLs), 7)
    mode(biosignMC) <- "character"
    
    checkIdentical(biosignMC, accMC)
    
  }
  
}


