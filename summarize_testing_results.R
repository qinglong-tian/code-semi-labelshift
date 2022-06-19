library(tidyverse)
read_in_data <- function(filename, dir) {
  rt <- str_match(filename, "ratio_(.*?)_")[2] %>% as.numeric()
  n <- str_match(filename, "_n(.*?)_ratio")[2] %>% as.numeric()
  beta1 <- str_match(filename, "beta1_(.*?)_.RDS")[2] %>% as.numeric()
  dat <- readRDS(paste(dir, filename, sep = ""))
  
  return(list(rt = rt, dat = dat, n = n, beta1 = beta1))
}
