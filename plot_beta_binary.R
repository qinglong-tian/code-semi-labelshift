source("summarize_beta_binary_results.R")

df_final %>% mutate(
  Method2 = factor(
    Method2,
    levels = c("Proposed", "Lipton"),
    labels = c("Proposed", "Lipton")
  ),
  Method = factor(
    Method,
    levels = c("logit", "rf", "svm", "mlp", "gbm", "nb"),
    labels = c("Logistic", "RF", "SVM", "MLP", "GBM", "NB")
  )
) -> df_final

df_final %>% filter(Method != "NB") -> df_final

df_final %>% filter(
  Type == "bias"
) %>% ggplot(
  aes(x = nVec, y = Value, col = Method2)
) + geom_line(
  aes(linetype = Method2)
) + facet_wrap(
  vars(Method)
) + geom_hline(
  yintercept = 0, linetype = 2
) + xlab(
  "n"
) + ylab(
  "Bias"
)

df_final %>% filter(
  Type == "mse"
) %>% ggplot(
  aes(x = nVec, y = Value, col = Method2)
) + geom_line(
  aes(linetype = Method2)
) + facet_wrap(
  vars(Method)
) + xlab(
  "n"
) + ylab(
  "MSE"
) + scale_y_continuous(
  trans='log10'
)
library(knitr)
collapse_rows_dt <- expand.grid(
  "n+m" = c('500', '1000', '1500', '2000', '2500', '3000'),
  "Model" = c("Logistic", "RF", "SVM", "MLP", "GBM")
)
collapse_rows_dt <- collapse_rows_dt[c("Model", "n+m")]

df_final %>% filter(
  Type == "bias", Method2 == "Proposed"
) %>% arrange(
  Method, nVec
) %>% select(Value) -> biasVec
collapse_rows_dt$Bias <- round(biasVec, digits = 3)

df_final %>% filter(
  Type == "sd", Method2 == "Proposed"
) %>% arrange(
  Method, nVec
) %>% select(Value) -> sdVec
collapse_rows_dt$Sd <- format(sdVec, digits = 3)

df_final %>% filter(
  Type == "se", Method2 == "Proposed"
) %>% arrange(
  Method, nVec
) %>% select(Value) -> seVec
collapse_rows_dt$Se <- format(seVec, digits = 3)

df_final %>% filter(
  Type == "cp", Method2 == "Proposed"
) %>% arrange(
  Method, nVec
) %>% select(Value) -> cpVec
collapse_rows_dt$CP <- format(cpVec, digits = 3)

library(kableExtra)

kbl(
  collapse_rows_dt,
  booktabs = T,
  align = 'c',
  format = 'latex'
) %>%
  collapse_rows(columns = 1:2, latex_hline = "major", valign = "middle")
