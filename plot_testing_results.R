library(tidyverse)
library(reshape2)
source("summarize_testing_results.R")

dat4plot %>% melt(
  id.vars = c("n", "ratio", "beta1"),
  value.name = "Rej",
  measure.vars = c("rejEff", "rejNaive"),
  variable.name = "Method"
) %>% mutate(
  Method = factor(
    Method,
    levels = c("rejEff", "rejNaive"),
    labels = c("Proposed", "Naive")
  ),
  Ratio = ifelse(
    ratio == 0.5,
    'm/n*"=0.5"',
    ifelse(ratio == 1, 'm/n*"=1"', 'm/n*"=1.5"')
  ),
  N = ifelse(n == 400,
             'n*"=400"',
             ifelse(n == 700, 'n*"=700"', 'n*"=1000"'))
) %>% mutate(N = factor(N, levels = c('n*"=400"', 'n*"=700"', 'n*"=1000"')),
             Ratio = factor(Ratio, levels = c('m/n*"=0.5"', 'm/n*"=1"', 'm/n*"=1.5"'))) %>% ggplot(aes(x = beta1, y = Rej, col = Method)) +
  geom_point(aes(shape = Method)) + geom_line(aes(linetype = Method)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  facet_grid(Ratio ~ N, labeller = label_parsed) + ylab("Rejection Rate")+xlab(expression(beta[1]))
