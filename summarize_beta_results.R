# Useful snippets
##############################
library(stringr)
source("data_generating_functions.R")

read_in_data <- function(filename, dir) {
  rt <- str_match(filename, "ratio_(.*?)_")[2] %>% as.numeric()
  n <- str_match(filename, "n_(.*?)_ratio")[2] %>% as.numeric()
  dat <- readRDS(paste(dir, filename, sep = ""))
  
  return(list(rt = rt, dat = dat, n = n))
}

colMedian <- function(v)
{
  apply(v, 2, median)
}

colSd <- function(v)
{
  apply(v, 2, sd)
}

colSquareMean <- function(v)
{
  apply(v, 2, function(x) {
    mean(x ^ 2)
  })
}

remove_outliers <- function(mat)
{
  IQRs <- apply(mat, 2, IQR)
  Bounds <- apply(mat, 2, quantile, prob = c(0.25, 0.75))
  Bounds[1, ] <- Bounds[1, ] - 1.5 * IQRs
  Bounds[2, ] <- Bounds[2, ] + 1.5 * IQRs
  
  mat_clean <- NULL
  removed <- NULL
  for (i in 1:nrow(mat))
  {
    bool_b1 <- (mat[i, 1] > Bounds[1, 1] & mat[i, 1] < Bounds[2, 1])
    bool_b2 <- (mat[i, 2] > Bounds[1, 2] & mat[i, 2] < Bounds[2, 2])
    
    if (bool_b1 & bool_b2)
    {
      mat_clean <- rbind(mat_clean, mat[i, ])
    }
    else
    {
      removed <- c(removed, i)
    }
  }
  
  return(list(mat_clean = mat_clean, removed = removed))
}

bias_mat <- function(mat, trueVal)
{
  mat - matrix(
    trueVal,
    ncol = length(trueVal),
    nrow = nrow(mat),
    byrow = T
  )
}
##############################
dir_to_dat <- "dat3/"
all_filenames <- list.files(dir_to_dat)

# Data Generating Parameters
Mu_YX <- c(2, 1, 1, 1)
SigMat_YX <- ar1_cor(4, 0.9)
SigMat_YX[-1,-1] <- ar1_cor(3, 0.3)
SigMat_YX <- SigMat_YX + 0.1 * diag(4)
SigMat_YX[1, 1] <- 1.44
Mu_Y_T <- 1.5
Sig_Y_T <- 1.5
Mu_Y_S <- Mu_YX[1]
Sig_Y_S <- sqrt(SigMat_YX[1, 1])

trueBetaRho <-
  Compute_Rho_Parameters(Mu_Y_T, Sig_Y_T, Mu_YX, SigMat_YX)
##############################
datBeta <- NULL
col_names <- c(
  "b1eff",
  "b2eff",
  "b1effbias",
  "b2effbias",
  "b1naive",
  "b2naive",
  "b1naivebias",
  "b2naivebias",
  "b1effse",
  "b2effse",
  "b1naivese",
  "b2naivese",
  "b1effmse",
  "b2effmse",
  "b1naivemse",
  "b2naivemse",
  "b1effpertsd",
  "b2effpertsd",
  "b1effformusd",
  "b2effformusd",
  "b1effformucp",
  "b2effformucp",
  "b1effpertcp",
  "b2effpertcp",
  "n",
  "ratio"
)
for (filename in all_filenames)
{
  infoList <- read_in_data(filename, dir_to_dat)
  rt <- infoList$rt
  n <- infoList$n
  result_list <- infoList$dat
  
  firstList <- result_list$Out1
  secondList <- result_list$Out2
  
  validVec <- sapply(firstList, function (x) {
    class(x) != "try-error"
  })
  
  first_list_clean <- vector("list", sum(validVec))
  second_list_clean <- vector("list", sum(validVec))
  
  index <- 1
  for (i in 1:length(validVec))
  {
    if (validVec[i])
    {
      first_list_clean[index] <- firstList[i]
      second_list_clean[index] <- secondList[i]
      index <- index + 1
    }
  }
  
  outVec <- NULL
  nameVec <- NULL
  
  # Beta Hat Efficient
  BetaHatEffMat <- t(sapply(first_list_clean, function (x) {
    x$BetaHatEff
  }))
  tmp <- remove_outliers(BetaHatEffMat)
  tmp$mat_clean -> BetaHatEffMat
  removed_eff <- tmp$removed
  
  outVec <- c(outVec, colMeans(BetaHatEffMat))

  # Bias Efficient
  
  outVec <- c(outVec, colMeans(BetaHatEffMat) - trueBetaRho)

  # Beta Hat Naive
  BetaHatNaiveMat <- t(sapply(first_list_clean, function (x) {
    x$BetaHatNaive
  }))
  remove_outliers(BetaHatNaiveMat)$mat_clean -> BetaHatNaiveMat
  
  outVec <- c(outVec, colMeans(BetaHatNaiveMat))

  # Bias Naive
  outVec <- c(outVec, colMeans(BetaHatNaiveMat) - trueBetaRho)

  # SE Eff
  outVec <- c(outVec, colSd(BetaHatEffMat))

  # SE Naive
  outVec <- c(outVec, colSd(BetaHatNaiveMat))

  # MSE Eff
  biasMatEff <- bias_mat(BetaHatEffMat, trueBetaRho)
  mseEff <- colSquareMean(biasMatEff)
  
  # MSE Naive
  biasMatNaive <- bias_mat(BetaHatNaiveMat, trueBetaRho)
  mseNaive <- colSquareMean(biasMatNaive)
  
  outVec <- c(outVec, mseEff, mseNaive)

  # Eff: Perturbation v.s., Formula
  
  PerturbSd <- t(sapply(second_list_clean, function (PertMat) {
    apply(PertMat, MARGIN = 2, sd)
  }))
  PerturbSd <- PerturbSd[-removed_eff,]
  outVec <- c(outVec, colMedian(PerturbSd))

  EffSdMat <- t(sapply(first_list_clean, function (x) {
    x$SdEff
  }))
  EffSdMat <- EffSdMat[-removed_eff,]
  outVec <- c(outVec, colMedian(EffSdMat))

  CPEffMat <- t(sapply(first_list_clean, function (x) {
    x$CPEff
  }))
  CPEffMat <- CPEffMat[-removed_eff,]
  outVec <- c(outVec, colMeans(CPEffMat))

  CPPert1 <- mean(
    (trueBetaRho[1] > BetaHatEffMat[, 1] - 1.96 * PerturbSd[, 1]) &
      (trueBetaRho[1] < BetaHatEffMat[, 1] + 1.96 * PerturbSd[, 1])
  )
  CPPert2 <- mean(
    (trueBetaRho[2] > BetaHatEffMat[, 2] - 1.96 * PerturbSd[, 2]) &
      (trueBetaRho[2] < BetaHatEffMat[, 2] + 1.96 * PerturbSd[, 2])
  )
  outVec <- c(outVec, CPPert1, CPPert2)

  outVec <- c(outVec, n, rt)
  outVec <- matrix(outVec, nrow = 1)
  
  datBeta <- rbind(datBeta, outVec)
}
colnames(datBeta) <- col_names
datBeta <- as.data.frame(datBeta)
