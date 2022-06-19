library(tidyverse)
read_in_data <- function(filename, dir) {
  rt <- str_match(filename, "ratio_(.*?)_")[2] %>% as.numeric()
  n <- str_match(filename, "_n(.*?)_ratio")[2] %>% as.numeric()
  beta1 <- str_match(filename, "beta1_(.*?)_.RDS")[2] %>% as.numeric()
  dat <- readRDS(paste(dir, filename, sep = ""))
  
  return(list(rt = rt, dat = dat, n = n, beta1 = beta1))
}

dir_to_dat <- "dat4/"
all_filenames <- list.files(dir_to_dat)

dat4plot <- matrix(nrow = length(all_filenames), ncol = 5)
i <- 1
for (fname in all_filenames)
{
  infoList <- read_in_data(fname, dir_to_dat)
  beta1 <- infoList$beta1
  n <- infoList$n
  ratio <- infoList$rt
  
  dat <- infoList$dat
  pEffVec <- sapply(dat, function(x) {x$pValEff})
  pNaiveVec <- sapply(dat, function(x) {x$pValNaive})
  
  rejEff <- mean(pEffVec <= 0.05, na.rm = T)
  rejNaive <- mean(pNaiveVec <= 0.05, na.rm = T)
  
  dat4plot[i,] <- c(n, ratio, beta1, rejEff, rejNaive)
  i <- i+1
}
colnames(dat4plot) <- c("n", "ratio", "beta1", "rejEff", "rejNaive")
dat4plot <- as.data.frame(dat4plot)
