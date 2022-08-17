library(parallel)
library(doParallel)
library(ranger)
source("data_generating_functions_binary.R")
source("estimation_functions_binary.R")

n <- 2000
B1 <- 200

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
# saveRDS(data_list_mc,
#         file = paste("binary_dat/MC_Data_List_n_", n, ".RDS", sep = ""))

# Model Selection

control <- trainControl(method = "cv",
                        number = 10,
                        search = "random")
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
  
  # Random Forest
  rfControl <- trainControl(method = "cv",
                            number = 10,
                            search = "grid")
  rftunegrid <-
    expand.grid(
      .mtry = c(1, 2),
      .min.node.size = c(100, 150, 200),
      .splitrule = c("extratrees", "gini")
    )
  rfFit2 <- train(
    Y ~ .,
    data = sData,
    method = "ranger",
    metric = "Accuracy",
    trControl = rfControl,
    tuneGrid = rftunegrid
  )
  
  dat.ranger <- ranger(
    Y ~ .,
    data = sData,
    probability = T,
    splitrule = rfFit2$finalModel$splitrule,
    mtry = rfFit2$finalModel$mtry,
    min.node.size = rfFit2$finalModel$min.node.size
  )
  rfProbs <-
    predict(dat.ranger, data = sData[, -1], type = "response")$predictions[, 2]
  rfProbt <-
    predict(dat.ranger, data = tData, type = "response")$predictions[, 2]
  
  # MLP
  mlpTrain <- train(
    Y ~ .,
    data = sData,
    method = "nnet",
    metric = metric,
    trControl = control
  )
  
  train(
    Y ~ .,
    data = sData,
    method = "nnet",
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
  
  nbProbs <- predict(nbFit, sData[, -1], type = "prob")[, 2]
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
  svmProbs <- predict(svmFit, sData[, -1], type = "prob")[, 2]
  svmProbt <- predict(svmFit, tData, type = "prob")[, 2]
  
  # GBM
  
  gbmControl <- trainControl(method = "cv",
                             number = 10,
                             search = "grid")
  gbmGrid <-
    expand.grid(
      .n.minobsinnode = c(3, 5, 10, 15, 20),
      .n.trees = c(800, 2000, 500, 1000, 1500),
      .interaction.depth = 1,
      .shrinkage = 0.01
    )
  
  train(
    Y ~ .,
    data = sData,
    method = "gbm",
    metric = "Accuracy",
    trControl = gbmControl,
    tuneGrid = gbmGrid
  ) -> gbmTrain
  
  gbmFit <- train(
    Y ~ .,
    data = sData,
    method = "gbm",
    tuneGrid = data.frame(gbmTrain$finalModel$tuneValue),
    trControl = trainControl(method = "none")
  )
  
  gbmProbs <-
    predict(gbmFit, as.matrix(sData[,-1]), type = "prob")[, 2]
  gbmProbt <-
    predict(gbmFit, as.matrix(tData), type = "prob")[, 2]
  
  return(
    list(
      logit_s = logitProbs,
      logit_t = logitProbt,
      rf_s = rfProbs,
      rf_t = rfProbt,
      mlp_s = mlpProbs,
      mlp_t = mlpProbt,
      nb_s = nbProbs,
      nb_t = nbProbt,
      svm_s = svmProbs,
      svm_t = svmProbt,
      xgb_s = gbmProbs,
      xgb_t = gbmProbt
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
  rf_s <- probTemp$rf_s
  rf_t <- probTemp$rf_t
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
  
  # Random Forest
  optim(
    par = trueBeta,
    fn = Compute_S_Eff_Sum,
    method = "BFGS",
    sData = sData,
    piVal = piVal,
    yFittedGivenX = rf_s,
    yFittedGivenX_tDat = rf_t
  ) -> optimOut
  c_ps <- E_S_RHO_Binary(optimOut$par, sData)
  
  Compute_SE_Binary(
    betaVal = optimOut$par,
    sData = sData,
    c_ps = c_ps,
    piVal = piVal,
    yFittedGivenX = rf_s,
    yFittedGivenX_tDat = rf_t
  ) -> se
  lwb <- optimOut$par - 1.96 * se
  upb <- optimOut$par + 1.96 * se
  cp <- trueBeta < upb & trueBeta > lwb
  
  beta_rf <- optimOut$par
  se_rf <- se
  cp_rf <- cp
  
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
    beta_rf,
    se_rf,
    cp_rf,
    beta_mlp,
    se_mlp,
    cp_mlp,
    beta_nb,
    se_nb,
    cp_nb,
    beta_svm,
    se_svm,
    cp_svm,
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
