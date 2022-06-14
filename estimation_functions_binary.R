# Functions to derive the efficient estimator
E_S_RHO_Binary <- function(betaVal, sData)
{
  yVec <- sData[,"Y"]
  prY1 <- mean(yVec)
  prY0 <- 1-prY1
  out <- prY0+prY1*exp(betaVal)
  
  return(out)
}


