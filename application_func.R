library(tidyverse)

Compute_Marginal_of_SOFA <- function(sData, tData)
{
  tSofa <- unique(tData$sofa)
  sSofa <- unique(sData$sofa)
  
  if (any(tSofa %in% sSofa == FALSE))
  {
    stop("The support of target data is not a subset of that of the source data.")
  }
  
  pSofaS <- table(sData$sofa) / length(sData$sofa)
  pSofat <- table(tData$sofa) / length(tData$sofa)
  sofaClassesS <- dimnames(pSofaS)[[1]]
  sofaClassesT <- dimnames(pSofat)[[1]]
  
  dfOut <- data.frame()
  for (i in 1:length(sofaClassesS))
  {
    sofaclass <- sofaClassesS[i]
    pt <- ifelse(sofaclass %in% sofaClassesT, pSofat[sofaclass], 0)
    ps <- pSofaS[sofaclass]
    ratio <- pt / ps
    dfOut[i, 1] <- sofaclass
    dfOut[i, 2] <- ratio
    dfOut[i, 3] <- ps
    dfOut[i, 4] <- pt
  }
  colnames(dfOut) <- c("SOFA", "Ratio", "Ps", "Pt")
  return(dfOut)
}

Compute_Pr_R1_Given_Y <- function(sData, tData)
{
  num_of_s <- nrow(sData)
  num_of_t <- nrow(tData)
  
  rMarginalRatio <- num_of_t / num_of_s
  df <- Compute_Marginal_of_SOFA(sData, tData)
  dff <- df %>% mutate(rProb = 1 / (Ratio * rMarginalRatio + 1))
  
  return(dff)
}

Compute_GLM_Pr_R1_Given_Y <- function(sData, tData, rSofaLogistic)
{
  dff <- Compute_Pr_R1_Given_Y(sData, tData)
  as.numeric(dff$SOFA) -> sofa_num
  predict(rSofaLogistic,
          newdata = data.frame(sofa = sofa_num),
          type = "response") -> predProb
  predProb <- unname(predProb)
  dff %>% mutate(rProbPred = predProb) -> dff
  
  return(dff)
}
