library(parallel)
library(doParallel)
library(ranger)
source("data_generating_functions_binary.R")
source("estimation_functions_binary.R")

n <- 2000
B1 <- 5

num_of_x <- 2
Mu_X <- rep(0, num_of_x)
Sigma_X <- matrix(nrow = num_of_x, ncol = num_of_x)
rho <- 0.8
for (i in 1:num_of_x)
{
  for (j in 1:num_of_x)
  {
    Sigma_X[i, j] <- rho ^ (abs(i - j))
  }
}
alphaVec <- c(0, rep(1, length(Mu_X)))
trueBeta <- -2
gammaVec <- c(0, -trueBeta)

dir.create(file.path("binary_dat/"))

###############
# Preparation
###############

# Data

mclapply(1:B1, function(x)
{
  Generate_Binary_Data(n, Mu_X, Sigma_X, alphaVec, gammaVec)
},
mc.cores = detectCores()) -> data_list_mc

saveRDS(data_list_mc,
        file = paste("binary_dat/MC_Data_List_n_", n, ".RDS", sep = ""))

# Model Selection

control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 3,
  search = "random"
)
metric <- "Accuracy"

# Estimating Probability

lapply(data_list_mc, function(dat)
{
  tData <- as.data.frame(dat$tData)
  sData <- as.data.frame(dat$sData)
  sData[, "Y"] <-
    factor(sData[, "Y"],
           levels = c(0, 1),
           labels = c("no", "yes"))
  
  # Logistic
  fittedVal <- Generate_Y_Given_X(sData[,-1], tData, sData)
  logitProbs <- fittedVal$prob_s
  logitProbt <- fittedVal$prob_t
  
  # MLP
  
  mlpTrain <- train(
    Y ~ .,
    data = sData,
    method = "mlp",
    metric = metric,
    trControl = control
  )
  
  train(
    Y ~ .,
    data = sData,
    method = "mlp",
    tuneGrid = mlpTrain$finalModel$tuneValue,
    trControl = trainControl(method = "none"),
  ) -> mlpFit
  
  mlpProbs <-
    predict(mlpFit, newdata = sData[, -1], type = "prob")[, 2]
  mlpProbt <- predict(mlpFit, newdata = tData, type = "prob")[, 2]
  
  # Naive Bayes
  
  nbFit <- train(
    Y ~ .,
    data = sData,
    method = 'nb',
    metric = metric,
    trControl = control
  )
  
  nbProbs <- predict(nbFit, sData[,-1], type = "prob")[, 2]
  nbProbt <- predict(nbFit, tData, type = "prob")[, 2]
  
  # SVM
  
  svmFit <- train(
    Y ~ .,
    data = sData,
    method = 'svmLinear2',
    metric = metric,
    trControl = control,
    probability = TRUE
  )
  svmProbs <- predict(svmFit, sData[,-1], type = "prob")[, 2]
  svmProbt <- predict(svmFit, tData, type = "prob")[, 2]
  
  # XGBoost
  
  train(
    Y ~ .,
    data = sData,
    method = "xgbTree",
    metric = metric,
    trControl = control
  ) -> xgbFit
  
  xgboostFit <- train(
    Y ~ .,
    data = sData,
    method = "xgbTree",
    tuneGrid = data.frame(xgbFit$finalModel$tuneValue),
    trControl = trainControl(method = "none")
  )
  
  xgbProbs <-
    predict(xgboostFit, as.matrix(sData[,-1]), type = "prob")[, 2]
  xgbProbt <-
    predict(xgboostFit, as.matrix(tData), type = "prob")[, 2]
  
  return(
    list(
      logit_s = logitProbs,
      logit_t = logitProbt,
      mlp_s = mlpProbs,
      mlp_t = mlpProbt,
      nb_s = nbProbs,
      nb_t = nbProbt,
      svm_s = svmProbs,
      svm_t = svmProbt,
      xgb_s = xgbProbs,
      xgb_t = xgbProbt
    )
  )
}) -> probList

##########################
# Estimation & Inference
##########################

estimation_inference <- mclapply(1:B1, function(i)
{
  YXR <- data_list_mc[[i]]
  tData <- YXR$tData
  sData <- YXR$sData
  sData <- as.data.frame(sData)
  piVal <- nrow(sData) / n
  
  probTemp <- probList[[i]]
  logit_s <- probTemp$logit_s
  logit_t <- probTemp$logit_t
  mlp_s <- probTemp$mlp_s
  mlp_t <- probTemp$mlp_t
  nb_s <- probTemp$nb_s
  nb_t <- probTemp$nb_t
  svm_s <- probTemp$svm_s
  svm_t <- probTemp$svm_t
  xgb_s <- probTemp$xgb_s
  xgb_t <- probTemp$xgb_t
  
  # Logit
  optim(
    par = trueBeta,
    fn = Compute_S_Eff_Sum,
    method = "BFGS",
    sData = sData,
    piVal = piVal,
    yFittedGivenX = logit_s,
    yFittedGivenX_tDat = logit_t
  ) -> optimOut
  c_ps <- E_S_RHO_Binary(optimOut$par, sData)
  
  Compute_SE_Binary(
    betaVal = optimOut$par,
    sData = sData,
    c_ps = c_ps,
    piVal = piVal,
    yFittedGivenX = logit_s,
    yFittedGivenX_tDat = logit_t
  ) -> se
  lwb <- optimOut$par - 1.96 * se
  upb <- optimOut$par + 1.96 * se
  cp <- trueBeta < upb & trueBeta > lwb
  
  beta_logit <- optimOut$par
  se_logit <- se
  cp_logit <- cp
  
  # MLP
  optim(
    par = -1.5,
    fn = Compute_S_Eff_Sum,
    method = "BFGS",
    sData = sData,
    piVal = piVal,
    yFittedGivenX = mlp_s,
    yFittedGivenX_tDat = mlp_t
  ) -> optimOut
  c_ps <- E_S_RHO_Binary(optimOut$par, sData)
  
  Compute_SE_Binary(
    betaVal = optimOut$par,
    sData = sData,
    c_ps = c_ps,
    piVal = piVal,
    yFittedGivenX = mlp_s,
    yFittedGivenX_tDat = mlp_t
  ) -> se
  lwb <- optimOut$par - 1.96 * se
  upb <- optimOut$par + 1.96 * se
  cp <- trueBeta < upb & trueBeta > lwb
  
  beta_mlp <- optimOut$par
  se_mlp <- se
  cp_mlp <- cp
  
  # Naive Bayes
  optim(
    par = -1.5,
    fn = Compute_S_Eff_Sum,
    method = "BFGS",
    sData = sData,
    piVal = piVal,
    yFittedGivenX = nb_s,
    yFittedGivenX_tDat = nb_t
  ) -> optimOut
  c_ps <- E_S_RHO_Binary(optimOut$par, sData)
  
  Compute_SE_Binary(
    betaVal = optimOut$par,
    sData = sData,
    c_ps = c_ps,
    piVal = piVal,
    yFittedGivenX = nb_s,
    yFittedGivenX_tDat = nb_t
  ) -> se
  lwb <- optimOut$par - 1.96 * se
  upb <- optimOut$par + 1.96 * se
  cp <- trueBeta < upb & trueBeta > lwb
  
  beta_nb <- optimOut$par
  se_nb <- se
  cp_nb <- cp
  
  # SVM
  optim(
    par = -1.5,
    fn = Compute_S_Eff_Sum,
    method = "BFGS",
    sData = sData,
    piVal = piVal,
    yFittedGivenX = svm_s,
    yFittedGivenX_tDat = svm_t
  ) -> optimOut
  c_ps <- E_S_RHO_Binary(optimOut$par, sData)
  
  Compute_SE_Binary(
    betaVal = optimOut$par,
    sData = sData,
    c_ps = c_ps,
    piVal = piVal,
    yFittedGivenX = svm_s,
    yFittedGivenX_tDat = svm_t
  ) -> se
  lwb <- optimOut$par - 1.96 * se
  upb <- optimOut$par + 1.96 * se
  cp <- trueBeta < upb & trueBeta > lwb
  
  beta_svm <- optimOut$par
  se_svm <- se
  cp_svm <- cp
  
  # xgBoost
  optim(
    par = -1.5,
    fn = Compute_S_Eff_Sum,
    method = "BFGS",
    sData = sData,
    piVal = piVal,
    yFittedGivenX = xgb_s,
    yFittedGivenX_tDat = xgb_t
  ) -> optimOut
  c_ps <- E_S_RHO_Binary(optimOut$par, sData)
  
  Compute_SE_Binary(
    betaVal = optimOut$par,
    sData = sData,
    c_ps = c_ps,
    piVal = piVal,
    yFittedGivenX = xgb_s,
    yFittedGivenX_tDat = xgb_t
  ) -> se
  lwb <- optimOut$par - 1.96 * se
  upb <- optimOut$par + 1.96 * se
  cp <- trueBeta < upb & trueBeta > lwb
  
  beta_bl <- optimOut$par
  se_bl <- se
  cp_bl <- cp
  
  ValVec <- c(
    beta_logit,
    se_logit,
    cp_logit,
    beta_mlp,
    se_mlp,
    cp_mlp,
    beta_nb,
    se_nb,
    cp_nb,
    beta_svm,
    se_svm,
    cp_svm
    beta_bl,
    se_bl,
    cp_bl
  )
  
  return(ValVec)
},
mc.cores = detectCores())

saveRDS(
  estimation_inference,
  file = paste("binary_dat/Estimation_And_Inference_n_", n, ".RDS", sep = "")
)

Sys.time() - t0
