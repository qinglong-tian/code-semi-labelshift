library(tidyverse)
trueBeta <- -2
remove_outliers_binary <- function(vec)
{
  one_iqr <- IQR(vec, na.rm = T)
  lwrq <- quantile(vec, 0.25, na.rm = T) - 2 * one_iqr
  uprq <- quantile(vec, 0.75, na.rm = T) + 2 * one_iqr
  
  removed <- is.nan(vec) | is.na(vec) | (vec < lwrq) | (vec > uprq)
  out <- vec[!removed]
  
  return(out)
}

dat_dir <- "binary_dat/"
files <- list.files(dat_dir)

TypeVec <- rep(c("beta", "se", "cp"), 6)
MethodVec <- rep(c("logit", "rf", "mlp", "nb", "svm", "gbm"), each = 3)

df_final <- NULL
for (file in files)
{
  tmp <- readRDS(paste(dat_dir, file, sep = ""))
  n <- str_match(file, "Inference_n_(.*?).RDS")[2] %>% as.numeric()
  tmpMat <- t(sapply(tmp, function(x) {x}))
  nVec <- n
  
  betaMat <- tmpMat[, c(1, 4, 7, 10, 13, 16)]
  betaVec <- apply(betaMat, MARGIN = 2, function(x)
  {
    xx <- remove_outliers_binary(x)
    mean(xx)
  })
  
  biasVec <- apply(betaMat, MARGIN = 2, function(x)
  {
    xx <- remove_outliers_binary(x)
    mean(xx) - trueBeta
  })
  
  sdVec <- apply(betaMat, MARGIN = 2, function(x)
  {
    xx <- remove_outliers_binary(x)
    sd(xx)
  })
  
  cpMat <- tmpMat[, c(3, 6, 9, 12, 15, 18)]
  cpVec <- colMeans(cpMat)
  
  seMat <- tmpMat[, c(2, 5, 8, 11, 14, 17)]
  seVec <- apply(seMat, MARGIN = 2, function(x)
  {
    xx <- remove_outliers_binary(x)
    mean(xx)
  })
  
  apply(betaMat, 2, function(x)
  {
    xx <- remove_outliers_binary(x)
    mean((xx-trueBeta)^2)
  }) -> mseVec
  
  ValVec <- c(c(t(cbind(betaVec, seVec, cpVec))), mseVec, biasVec, sdVec)

  df <- data.frame(
    Value = ValVec,
    Type = c(TypeVec, rep("mse", 6), rep("bias", 6), rep("sd", 6)),
    Method = c(MethodVec, rep(c("logit", "rf", "mlp", "nb", "svm", "gbm"), 3)),
    nVec = nVec
  )
  
  df_final <- rbind(df_final, df)
}

df_final$Method2 <- "Proposed"

dat_dir1 <- "binary_dat1/"
files1 <- list.files(dat_dir1)
for (file in files1)
{
  dfTmp <- NULL
  tmp <- readRDS(paste(dat_dir1, file, sep = ""))
  n <- str_match(file, "Lipton_n_(.*?).RDS")[2] %>% as.numeric()
  tmp <- t(sapply(tmp, function(x) {
    unlist(x)
  }))
  
  # Bias
  biasVec <- colMeans(tmp) - trueBeta
  # MSE
  trueBetaMat <-
    matrix(trueBeta, ncol = ncol(tmp), nrow = nrow(tmp))
  mseMat <- tmp - trueBetaMat
  mseVec <- colMeans(mseMat ^ 2)
  
  valueVec <- c(biasVec, mseVec)
  typeVec <- c(rep("bias", 6), rep("mse", 6))
  nVec <- rep(n, 12)
  methodVec <- rep(c("logit", "rf", "mlp", "nb", "svm", "gbm"), 2)
  method2Vec <- rep("Lipton", 12)
  
  dfTmp$Value <- valueVec
  dfTmp$Type <- typeVec
  dfTmp$Method <- methodVec
  dfTmp$nVec <- nVec
  dfTmp$Method2 <- method2Vec
  
  df_final <- rbind(df_final, dfTmp)
}
