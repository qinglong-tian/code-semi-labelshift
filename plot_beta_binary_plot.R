source("summarize_beta_binary_results.R")
library(ggpubr)

# scales::hue_pal()(5)

df_final %>% filter(
  Type == "beta"
)  %>% mutate(
  Bias = Value-trueBeta
)  %>% ggplot(
  aes(x = nVec, y = Bias, col = Method)
) + geom_point(
  aes(shape = Method)
) + geom_line(
  aes(linetype = Method)
) + xlab(
  "n"
)

df_final %>% filter(
  Type == "mse"
) %>% ggplot(
  aes(x= nVec, y = Value, col = Method)
) + geom_point(
  aes(shape = Method)
) + geom_line(
  aes(linetype = Method)
) + ylab(
  "MSE"
) + xlab(
  "n"
) -> fig_mse

df_final %>% filter(
  Type == "cp"
) %>% ggplot(
  aes(x = nVec, y = Value, col = Method)
) + geom_point(
  aes(shape = Method)
) + geom_line(
  aes(linetype = Method)
) + ylab(
  "Coverage Probability"
) + geom_hline(
  yintercept = 0.95, linetype = 2
) + xlab(
  "n"
) -> fig_cp

ggarrange(
  fig_bias,
  fig_mse,
  fig_cp,
  ncol = 1,
  nrow = 3
) %>% annotate_figure(top = "bl: boosted logistic regression; nb: nayes Bayes; rf: random forest")
