library(parallel)
library(Rcpp)
library(fastGHQuad)
source("data_generating_functions.R")
sourceCpp("fast_estimation_functions.cpp")
#############################################
# Parameters
B1 <- 2000

beta_2 <- 0

Mu_YX <- c(2, 1, 1, 1)
SigMat_YX <- ar1_cor(4, 0.9)
Mu_Y <- Mu_YX[1]
Var_Y <- SigMat_YX[1, 1]

X_Given_Y_List <- Compute_X_Given_Y(Mu_YX, SigMat_YX)
yx_dist <- Compute_Y_Given_X(Mu_YX, SigMat_YX)
coef_y_x_s_true <- c(yx_dist$Beta0, yx_dist$Beta1)
var_y_x_s_true <- yx_dist$VarYX
sigma_y_x_s_true <- sqrt(var_y_x_s_true)

gh_num <- 10
ghxw <- gaussHermiteData(gh_num)
xList <- ghxw$x
wList <- ghxw$w
#############################################
n_vector <- c(500, 1000, 1500)
mnratio_vector <- c(0.5, 1, 1.5)
beta_1_vector <- seq(from = -1,
                     to = 1,
                     length.out = 7)
#############################################
t0 <- Sys.time()
for (n in n_vector)
{
  for (mnratio in mnratio_vector)
  {
    for (beta_1 in beta_1_vector)
    {
      m <- mnratio * n
      trueBetaVal <- c(beta_1, beta_2)
      data_list_mc <- mclapply(1:B1, function(x) {
        Testing_Data_Gen(n, m, Mu_YX, SigMat_YX, trueBetaVal)
      }, mc.cores = detectCores())
      
      results_output <- mclapply(data_list_mc, function(dataList) {
        sData <- dataList$sDat
        tData <- dataList$tDat
        piVal <- n / (n + m)
        ispar <- T
        
        fityx <- lm(Y ~ ., data = as.data.frame(sData))
        coef_y_x_s_hat <- coef(fityx)
        sigma_y_x_s_hat <- sigma(fityx)
        
        Mu_Y_S_hat <- mean(sData[, "Y"])
        Sig_Y_S_hat <- sd(sData[, "Y"])
        
        parameters <-
          list(
            y_vec = sData[, "Y"],
            mu = Mu_Y_S_hat,
            sigma = Sig_Y_S_hat,
            xList = xList,
            wList = wList
          )
        
        betaHat <- optim(
          trueBetaVal,
          EstimateBetaFunc_CPP,
          sData = sData,
          tData = tData,
          piVal = piVal,
          tDat_ext = tData,
          coef_y_x_s = coef_y_x_s_hat,
          sigma_y_x_s = sigma_y_x_s_hat,
          ispar = ispar,
          parameters = parameters,
          xList = xList,
          wList = wList
        )$par
        
        HatCovMat <-
          EstimateBetaCovMat_CPP(
            beta_rho = betaHat,
            sData = sData,
            tData = tData,
            piVal = piVal,
            tDat_ext = tData,
            coef_y_x_s = coef_y_x_s,
            sigma_y_x_s = sigma_y_x_s,
            ispar = ispar,
            parameters = parameters,
            xList = xList,
            wList = wList
          )
        
        betaHatNaive <- tryCatch(
          optim(
            trueBetaVal,
            Estimate_Naive_Beta,
            sData = sData,
            tData = tData,
            piVal = piVal,
            ispar = ispar,
            parameters = parameters
          )$par
          ,
          error = function(e) {
            return(c(NA, NA))
          }
        )
        
        # Testing Procedure
        beta_H0 <- matrix(c(0, 0), ncol = 1)
        R <- diag(2)
        
        ## Efficient
        betaHat <- matrix(betaHat, ncol = 1)
        cshiq2 <-
          t(R %*% betaHat - beta_H0) %*% solve(HatCovMat / (n + m)) %*% (R %*% betaHat - beta_H0)
        pVal <- 1 - pchisq(cshiq2, df = 2)
        
        ## Naive
        
        if (is.na(betaHatNaive[1]))
        {
          pVal_naive <- NA
        }
        else
        {
          pVal_naive <- tryCatch({
            DerivMat <-
              Compute_Derivative_Mat(betaHatNaive, sData, tData, ispar, parameters)
            EMM <-
              Compute_Var_m(betaHatNaive,
                            sData,
                            tData,
                            ispar,
                            parameters,
                            piVal)
            HatCovMatNaive <-
              solve(DerivMat) %*% EMM %*% t(solve(DerivMat)) / (n + m)
            betaHatNaive <- matrix(betaHatNaive, ncol = 1)
            cshiq2_naive <-
              t(R %*% betaHatNaive - beta_H0) %*% solve(HatCovMatNaive) %*% (R %*% betaHatNaive - beta_H0)
            1 - pchisq(cshiq2_naive, df = 2)
          },
          error = function(e) {
            return(NA)
          })
        }
        
        return(
          list(
            BetaEff = betaHat,
            BetaNaive = betaHatNaive,
            CovMatEff = HatCovMat / (n + m),
            CovMatNaive = HatCovMatNaive,
            pValEff = pVal,
            pValNaive = pVal_naive
          )
        )
      }, mc.cores = detectCores())
      saveRDS(
        results_output,
        file = paste(
          "dat4/testing_n",
          n,
          "_ratio_",
          mnratio,
          "_beta1",
          beta_1,
          "_.RDS",
          sep = ""
        )
      )
    }
  }
}
#############################################
Sys.time() - t0