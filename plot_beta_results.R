library(tidyverse)
library(reshape2)
library(ggpubr)
library(kableExtra)

# Load the results
source("summarize_beta_results.R")

# Plot the Coverage
datBeta %>%
  select(b1effformucp, b2effformucp, b1effpertcp, b2effpertcp, n, ratio) %>%
  melt(
    value.name = "Coverage",
    measure.vars = c("b1effformucp", "b2effformucp", "b1effpertcp", "b2effpertcp"),
    id.vars = c("n", "ratio")
  ) %>% mutate(parameter = ifelse(
    variable %in% c("b1effformucp", "b1effpertcp"),
    "beta[1]",
    "beta[2]"
  )) %>% mutate(Method = ifelse(
    variable %in% c("b1effformucp", "b2effformucp"),
    "Formula",
    "Perturbation"
  )) %>% mutate(Ratio = ifelse(
    ratio == 0.5,
    'm/n*"=0.5"',
    ifelse(ratio == 1, 'm/n*"=1"', 'm/n*"=1.5"')
  )) %>% mutate(Ratio = factor(Ratio, levels = c('m/n*"=0.5"', 'm/n*"=1"', 'm/n*"=1.5"'))) %>%
  ggplot(aes(x = n, y = Coverage, col = Method)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(shape = Method)) +
  facet_grid(parameter ~ Ratio, labeller = label_parsed) +
  geom_hline(yintercept = 0.95, linetype = "dashed") + ylim(c(0.9, 1)) +
  ylab("Coverage\ Probability")

# Plot the Bias
datBeta %>% select(b1effbias, b2effbias, b1naivebias, b2naivebias, n, ratio) %>% melt(
  value.name = "Bias",
  meausre.vars = c("b1effbias", "b2effbias", "b1naivebias", "b2naivebias"),
  id.vars = c("n", "ratio")
) %>% mutate(parameter = ifelse(
  variable %in% c("b1effbias", "b1naivebias"),
  "beta[1]",
  "beta[2]"
)) %>% mutate(Method = ifelse(variable %in% c("b1effbias", "b2effbias"), "Proposed", "Naive")) %>% mutate(Ratio = ifelse(
  ratio == 0.5,
  'm/n*"=0.5"',
  ifelse(ratio == 1, 'm/n*"=1"', 'm/n*"=1.5"')
)) %>% mutate(Ratio = factor(Ratio, levels = c('m/n*"=0.5"', 'm/n*"=1"', 'm/n*"=1.5"')),
              Method = factor(Method, levels = c("Proposed", "Naive"))) %>% mutate(Method = as.factor(Method)) %>%
  ggplot(aes(x = n, y = Bias, col = Method)) + geom_line(aes(linetype = Method)) +
  geom_point(aes(shape = Method)) +  facet_grid(parameter ~ Ratio, labeller = label_parsed, scales =
                                                  "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed") -> biasfig1

# MSE Ratio
datBeta %>% mutate(b1MSERatio = b1effmse / b1naivemse,
                   b2MSERatio = b2effmse / b2naivemse) %>% select(b1MSERatio, b2MSERatio, n, ratio) %>% melt(
                     value.name = "MSE",
                     measure.vars = c("b1MSERatio", "b2MSERatio"),
                     id.vars = c("n", "ratio")
                   ) %>% mutate(
                     parameter = ifelse(variable == "b1MSERatio", "beta[1]", "beta[2]"),
                     Ratio = ifelse(
                       ratio == 0.5,
                       'm/n*"=0.5"',
                       ifelse(ratio == 1, 'm/n*"=1"', 'm/n*"=1.5"')
                     )
                   )%>% mutate(Ratio = factor(Ratio, levels = c('m/n*"=0.5"', 'm/n*"=1"', 'm/n*"=1.5"'))) %>% ggplot(aes(x = n, y = MSE)) + geom_line() + geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed") + facet_grid(parameter ~ Ratio, labeller = label_parsed) +
  ylab("MSE Ratio") + ylim(c(0, 1))

# MSE
datBeta %>% select(b1effmse, b1naivemse, b2effmse, b2naivemse, n, ratio) %>% melt(
  value.name = "MSE",
  measure.vars = c("b1effmse", "b1naivemse", "b2effmse", "b2naivemse"),
  id.vars = c("n", "ratio")
) %>% mutate(
  parameter = ifelse(variable %in% c("b1effmse", "b1naivemse"), "beta[1]", "beta[2]"),
  Ratio = ifelse(
    ratio == 0.5,
    'm/n*"=0.5"',
    ifelse(ratio == 1, 'm/n*"=1"', 'm/n*"=1.5"')
  ),
  Method = ifelse(variable %in% c("b1effmse", "b2effmse"), "Proposed", "Naive")
) %>% mutate(Ratio = factor(Ratio, levels = c('m/n*"=0.5"', 'm/n*"=1"', 'm/n*"=1.5"')),
             Method = factor(Method, levels = c("Proposed", "Naive"))) %>% ggplot(aes(x = n, y = MSE, col = Method)) +
  geom_line(aes(linetype = Method)) + geom_point(aes(shape = Method)) + facet_grid(parameter ~ Ratio, labeller = label_parsed, scales =
                                                                           "free_y") -> msefig1

ggarrange(biasfig1,
          msefig1,
          ncol = 1,
          labels = c("(a)", "(b)"))

datBeta %>% select(
  ratio,
  n,
  b1effbias,
  b1effse,
  b1effpertsd,
  b1effpertcp,
  b2effbias,
  b2effse,
  b2effpertsd,
  b2effpertcp
) %>% arrange(ratio, n) -> latexTable
latexTable <- latexTable[,-1]

latexTable <- cbind(latexTable, latexTable2[,-c(1,2)])
kbl(latexTable,
    booktabs = T,
    digits = 3,
    format = "latex") %>% kable_styling(position = "center") %>%
  pack_rows("m/n=0.5", 1, 10) %>%
  pack_rows("m/n=1", 11, 20) %>%
  pack_rows("m/n=1.5", 21, 30) %>%
  add_header_above(c(
    " ",
    "beta_1" = 4,
    "beta_2" = 4,
    "btheta" = 4
  ))

# SE and SD
datBeta %>% select(
  b1effse,
  b2effse,
  b1naivese,
  b2naivese,
  b1effpertsd,
  b2effpertsd,
  b1effformusd,
  b2effformusd,
  n,
  ratio
) %>% melt(
  value.name = "SE",
  measure.vars = c(
    "b1effse",
    "b2effse",
    "b1naivese",
    "b2naivese",
    "b1effpertsd",
    "b2effpertsd",
    "b1effformusd",
    "b2effformusd"
  ),
  id.vars = c("n", "ratio")
) %>% mutate(
  parameter = ifelse(
    variable %in% c("b1effse", "b1naivese", "b1effpertsd", "b1effformusd"),
    "beta[1]",
    "beta[2]"
  ),
  Ratio = ifelse(
    ratio == 0.5,
    'm/n*"=0.5"',
    ifelse(ratio == 1, 'm/n*"=1"', 'm/n*"=1.5"')
  ),
  Type = ifelse(
    variable %in% c("b1effse", "b2effse"),
    "Proposed (Empirical)",
    ifelse(
      variable %in% c("b1naivese", "b2naivese"),
      "Naive (Empirical)",
      ifelse(
        variable %in% c("b1effpertsd", "b2effpertsd"),
        "Eff. Perturbation (Est.)",
        "Eff. Formula (Est.)"
      )
    )
  )
) %>% mutate(
  Ratio = factor(Ratio, levels = c('m/n*"=0.5"', 'm/n*"=1"', 'm/n*"=1.5"')),
  Type = factor(Type, levels = c("Proposed (Empirical)", "Naive (Empirical)", "Eff. Perturbation (Est.)", "Eff. Formula (Est.)"))
) %>% ggplot(
  aes(x = n, y = SE, col = Type)
) + geom_line(
  aes(linetype = Type)
) + geom_point(
  aes(shape = Type)
) + facet_grid(
  parameter ~ Ratio,
  labeller = label_parsed,
  scales="free_y"
) + ylab("Standard Error")

# SE
datBeta %>% select(
  b1effse,
  b2effse,
  b1naivese,
  b2naivese,
  b1effpertsd,
  b2effpertsd,
  b1effformusd,
  b2effformusd,
  n,
  ratio
) %>% melt(
  value.name = "SE",
  measure.vars = c(
    "b1effse",
    "b2effse",
    "b1naivese",
    "b2naivese",
    "b1effpertsd",
    "b2effpertsd",
    "b1effformusd",
    "b2effformusd"
  ),
  id.vars = c("n", "ratio")
) %>% mutate(
  parameter = ifelse(
    variable %in% c("b1effse", "b1naivese", "b1effpertsd", "b1effformusd"),
    "beta[1]",
    "beta[2]"
  ),
  Ratio = ifelse(
    ratio == 0.5,
    'm/n*"=0.5"',
    ifelse(ratio == 1, 'm/n*"=1"', 'm/n*"=1.5"')
  ),
  Type = ifelse(
    variable %in% c("b1effse", "b2effse"),
    "Proposed",
    ifelse(
      variable %in% c("b1naivese", "b2naivese"),
      "Naive",
      ifelse(
        variable %in% c("b1effpertsd", "b2effpertsd"),
        "Eff. Perturbation (Est.)",
        "Eff. Formula (Est.)"
      )
    )
  )
) %>% mutate(
  Ratio = factor(Ratio, levels = c('m/n*"=0.5"', 'm/n*"=1"', 'm/n*"=1.5"')),
  Type = factor(Type, levels = c("Naive", "Proposed", "Eff. Perturbation (Est.)", "Eff. Formula (Est.)"))
) -> tmp

tmp %>% filter(Type %in% c("Proposed", "Naive")) %>%  ggplot(
  aes(x = n, y = SE, col = Type)
) + geom_line(
  aes(linetype = Type)
) + geom_point(
  aes(shape = Type)
) + facet_grid(
  parameter ~ Ratio,
  labeller = label_parsed,
  scales="free_y"
) + ylab("Standard Deviation")
