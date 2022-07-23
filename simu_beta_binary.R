library(parallel)
source("data_generating_functions_binary.R")
source("estimation_functions_binary.R")

n <- 3000
Mu_X <- c(1, -1, 2)
Sigma_X <- diag(3)
alphaVec <- c(-1, 1, -2, -1)
trueBeta <- -2
gammaVec <- c(-1, -trueBeta)
Method <- "logit"
B1 <- 500

results <- mclapply(1:B1, function(x)
{
  YXR <- Generate_Binary_Data(n, Mu_X, Sigma_X, alphaVec, gammaVec)
  tData <- YXR$tData
  sData <- YXR$sData
  piVal <- nrow(sData) / n
  
  fittedVal <- Generate_Y_Given_X(sData[, -1], tData, sData, Method)
  probVecs <- fittedVal$prob_s
  probVect <- fittedVal$prob_t
  
  optim(
    par = trueBeta,
    fn = Compute_S_Eff_Sum,
    method = "BFGS",
    sData = sData,
    piVal = piVal,
    yFittedGivenX = probVecs,
    yFittedGivenX_tDat = probVect
  ) -> optim_beta
  
  c_ps <- E_S_RHO_Binary(optim_beta$par, sData)
  
  Compute_SE_Binary(
    betaVal = optim_beta$par,
    sData = sData,
    c_ps = c_ps,
    piVal = piVal,
    yFittedGivenX = probVecs,
    yFittedGivenX_tDat = probVect
  ) -> se
  
  lb <- optim_beta$par - 1.96 * se
  ub <- optim_beta$par + 1.96 * se
  
  cp <- trueBeta <= ub & trueBeta >= lb
  
  return(list(optim_beta$par,
              se,
              cp))
},
mc.cores = detectCores())

mean(sapply(results, function(x) {
  x[[3]]
}))

###############

B2 <- 100

t0 <- Sys.time()

data_list_mc <- mclapply(1:B2, function(x)
{
  Generate_Binary_Data(n, Mu_X, Sigma_X, alphaVec, gammaVec)
},
mc.cores = detectCores()
)

library(doParallel)
cl <- makePSOCKcluster(16)
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
  
  probVecs <- predict(rfFit, sData[,-1], type = "prob")[,2]
  probVect <- predict(rfFit, tData, type = "prob")[,2]
  
  return(list(
    probs = probVecs,
    probt = probVect
  ))
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
  
  probVecs <- predict(rfFit, sData[,-1], type = "prob")[,2]
  probVect <- predict(rfFit, tData, type = "prob")[,2]
  
  return(list(
    probs = probVecs,
    probt = probVect
  ))
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
  
  probVecs <- predict(rfFit, sData[,-1], type = "prob")[,2]
  probVect <- predict(rfFit, tData, type = "prob")[,2]
  
  return(list(
    probs = probVecs,
    probt = probVect
  ))
}) -> lb_fit_list

stopCluster(cl)

Sys.time() - t0