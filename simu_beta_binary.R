library(parallel)
library(doParallel)
source("data_generating_functions_binary.R")
source("estimation_functions_binary.R")

n <- 1000

Mu_X <- c(1,-1, 2)
Sigma_X <- diag(3)
alphaVec <- c(-1, 1,-2,-1)
trueBeta <- -2
gammaVec <- c(-1,-trueBeta)
B1 <- 100

t0 <- Sys.time()

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

saveRDS(data_list_mc, file = paste("binary_dat/MC_Data_List_n_", n, ".RDS", sep = ""))

# Fitted Model

mclapply(data_list_mc, function(YXR)
{
  tData <- YXR$tData
  sData <- YXR$sData
  
  fittedVal <- Generate_Y_Given_X(sData[,-1], tData, sData)
  probVecs <- fittedVal$prob_s
  probVect <- fittedVal$prob_t
  
  return(list(probs = probVecs,
              probt = probVect))
},
mc.cores = detectCores()) -> logit_fit_list

mclapply(data_list_mc, function(YXR)
{
  tData <- YXR$tData
  sData <- YXR$sData
  
  fittedVal <-
    Generate_Y_Given_X(sData[,-1], tData, sData, "probit")
  probVecs <- fittedVal$prob_s
  probVect <- fittedVal$prob_t
  
  return(list(probs = probVecs,
              probt = probVect))
},
mc.cores = detectCores()) -> probit_fit_list

cl <- makePSOCKcluster(detectCores() - 4)
registerDoParallel(cl)

lapply(data_list_mc, function(YXR)
{
  control <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 3,
    search = "random"
  )
  
  metric <- "Accuracy"
  sData <- YXR$sData
  sData <- as.data.frame(sData)
  sData[, "Y"] <- as.factor(sData[, "Y"])
  tData <- YXR$tData
  
  rfFit <- train(
    Y ~ .,
    data = sData,
    method = "rf",
    metric = metric,
    trControl = control
  )
  
  probVecs <- predict(rfFit, sData[,-1], type = "prob")[, 2]
  probVect <- predict(rfFit, tData, type = "prob")[, 2]
  
  return(list(probs = probVecs,
              probt = probVect))
}) -> rf_fit_list

lapply(data_list_mc, function(YXR)
{
  control <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 3,
    search = "random"
  )
  
  metric <- "Accuracy"
  sData <- YXR$sData
  sData <- as.data.frame(sData)
  sData[, "Y"] <- as.factor(sData[, "Y"])
  tData <- YXR$tData
  rfFit <- train(
    Y ~ .,
    data = sData,
    method = "nb",
    metric = metric,
    trControl = control
  )
  
  probVecs <- predict(rfFit, sData[,-1], type = "prob")[, 2]
  probVect <- predict(rfFit, tData, type = "prob")[, 2]
  
  return(list(probs = probVecs,
              probt = probVect))
}) -> nb_fit_list

lapply(data_list_mc, function(YXR)
{
  control <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 3,
    search = "random"
  )
  
  metric <- "Accuracy"
  sData <- YXR$sData
  sData <- as.data.frame(sData)
  sData[, "Y"] <- as.factor(sData[, "Y"])
  tData <- YXR$tData
  rfFit <- train(
    Y ~ .,
    data = sData,
    method = "LogitBoost",
    metric = metric,
    trControl = control
  )
  
  probVecs <- predict(rfFit, sData[,-1], type = "prob")[, 2]
  probVect <- predict(rfFit, tData, type = "prob")[, 2]
  
  return(list(probs = probVecs,
              probt = probVect))
}) -> lb_fit_list

stopCluster(cl)

saveRDS(
  list(
    logit = logit_fit_list,
    probit = probit_fit_list,
    rf = rf_fit_list,
    nb = nb_fit_list,
    lb = lb_fit_list
  ),
  file = paste("binary_dat/Fitted_Models_n_", n, ".RDS", sep = "")
)

##########################
# Estimation & Inference
##########################

estimation_inference <- mclapply(1:B1, function(i)
{
  YXR <- data_list_mc[[i]]
  tData <- YXR$tData
  sData <- YXR$sData
  piVal <- nrow(sData) / n
  
  # Logit
  fittedVal <- logit_fit_list[[i]]
  probVecs <- fittedVal$probs
  probVect <- fittedVal$probt
  
  optim(
    par = trueBeta,
    fn = Compute_S_Eff_Sum,
    method = "BFGS",
    sData = sData,
    piVal = piVal,
    yFittedGivenX = probVecs,
    yFittedGivenX_tDat = probVect
  ) -> optimOut
  c_ps <- E_S_RHO_Binary(optimOut$par, sData)
  
  Compute_SE_Binary(
    betaVal = optimOut$par,
    sData = sData,
    c_ps = c_ps,
    piVal = piVal,
    yFittedGivenX = probVecs,
    yFittedGivenX_tDat = probVect
  ) -> se
  lwb <- optimOut$par - 1.96 * se
  upb <- optimOut$par + 1.96 * se
  cp <- trueBeta < upb & trueBeta > lwb
  
  beta_logit <- optimOut$par
  se_logit <- se
  cp_logit <- cp
  
  # Probit
  fittedVal <- probit_fit_list[[i]]
  probVecs <- fittedVal$probs
  probVect <- fittedVal$probt
  
  optim(
    par = trueBeta,
    fn = Compute_S_Eff_Sum,
    method = "BFGS",
    sData = sData,
    piVal = piVal,
    yFittedGivenX = probVecs,
    yFittedGivenX_tDat = probVect
  ) -> optimOut
  c_ps <- E_S_RHO_Binary(optimOut$par, sData)
  
  Compute_SE_Binary(
    betaVal = optimOut$par,
    sData = sData,
    c_ps = c_ps,
    piVal = piVal,
    yFittedGivenX = probVecs,
    yFittedGivenX_tDat = probVect
  ) -> se
  lwb <- optimOut$par - 1.96 * se
  upb <- optimOut$par + 1.96 * se
  cp <- trueBeta < upb & trueBeta > lwb
  
  beta_probit <- optimOut$par
  se_probit <- se
  cp_probit <- cp
  
  # Random Forest
  fittedVal <- rf_fit_list[[i]]
  probVecs <- fittedVal$probs
  probVect <- fittedVal$probt
  
  optim(
    par = trueBeta,
    fn = Compute_S_Eff_Sum,
    method = "BFGS",
    sData = sData,
    piVal = piVal,
    yFittedGivenX = probVecs,
    yFittedGivenX_tDat = probVect
  ) -> optimOut
  c_ps <- E_S_RHO_Binary(optimOut$par, sData)
  
  Compute_SE_Binary(
    betaVal = optimOut$par,
    sData = sData,
    c_ps = c_ps,
    piVal = piVal,
    yFittedGivenX = probVecs,
    yFittedGivenX_tDat = probVect
  ) -> se
  lwb <- optimOut$par - 1.96 * se
  upb <- optimOut$par + 1.96 * se
  cp <- trueBeta < upb & trueBeta > lwb
  
  beta_rf <- optimOut$par
  se_rf <- se
  cp_rf <- cp
  
  # Naive Bayes
  fittedVal <- nb_fit_list[[i]]
  probVecs <- fittedVal$probs
  probVect <- fittedVal$probt
  
  optim(
    par = trueBeta,
    fn = Compute_S_Eff_Sum,
    method = "BFGS",
    sData = sData,
    piVal = piVal,
    yFittedGivenX = probVecs,
    yFittedGivenX_tDat = probVect
  ) -> optimOut
  c_ps <- E_S_RHO_Binary(optimOut$par, sData)
  
  Compute_SE_Binary(
    betaVal = optimOut$par,
    sData = sData,
    c_ps = c_ps,
    piVal = piVal,
    yFittedGivenX = probVecs,
    yFittedGivenX_tDat = probVect
  ) -> se
  lwb <- optimOut$par - 1.96 * se
  upb <- optimOut$par + 1.96 * se
  cp <- trueBeta < upb & trueBeta > lwb
  
  beta_nb <- optimOut$par
  se_nb <- se
  cp_nb <- cp
  
  # Boosted Logistic Regression
  fittedVal <- lb_fit_list[[i]]
  probVecs <- fittedVal$probs
  probVect <- fittedVal$probt
  
  optim(
    par = trueBeta,
    fn = Compute_S_Eff_Sum,
    method = "BFGS",
    sData = sData,
    piVal = piVal,
    yFittedGivenX = probVecs,
    yFittedGivenX_tDat = probVect
  ) -> optimOut
  c_ps <- E_S_RHO_Binary(optimOut$par, sData)
  
  Compute_SE_Binary(
    betaVal = optimOut$par,
    sData = sData,
    c_ps = c_ps,
    piVal = piVal,
    yFittedGivenX = probVecs,
    yFittedGivenX_tDat = probVect
  ) -> se
  lwb <- optimOut$par - 1.96 * se
  upb <- optimOut$par + 1.96 * se
  cp <- trueBeta < upb & trueBeta > lwb
  
  beta_bl <- optimOut$par
  se_bl <- se
  cp_bl <- cp
  
  #
  ValVec <- c(
    beta_logit,
    se_logit,
    cp_logit,
    beta_probit,
    se_probit,
    cp_probit,
    beta_rf,
    se_rf,
    cp_rf,
    beta_nb,
    se_nb,
    cp_nb,
    beta_bl,
    se_bl,
    cp_bl
  )
  
  return(ValVec)
},
mc.cores = detectCores())

saveRDS(estimation_inference,
        file = paste("binary_dat/Estimation_And_Inference_n_", n, ".RDS", sep = ""))

ValVec <-
  rowMeans(sapply(estimation_inference, function(m) {
    m[!is.finite(m)] <- NA
    return(m)
  }), na.rm = T)
TypeVec <- rep(c("beta", "se", "cp"), 5)
MethodVec <- rep(c("logit", "probit", "rf", "nb", "bl"), each = 3)
nVec <- rep(n, 15)

df <- data.frame(
  Value = ValVec,
  Type = TypeVec,
  Method = MethodVec,
  nVec = nVec
)
saveRDS(df, file = paste("binary_dat/df_n_", n, ".RDS", sep = ""))
Sys.time() - t0
