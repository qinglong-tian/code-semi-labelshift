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
  
  ValVec <- c(c(t(cbind(betaVec, seVec, cpVec))), mseVec)

  df <- data.frame(
    Value = ValVec,
    Type = c(TypeVec, rep("mse", 6)),
    Method = c(MethodVec, c("logit", "rf", "mlp", "nb", "svm", "gbm")),
    nVec = nVec
  )
  
  df_final <- rbind(df_final, df)
}

sapply(lipton_method_inference, function(x) {unlist(x)}) %>% rowMeans()
