library(Rcpp)
library(brms)
library(readxl)
library(rstan)
library(tidybayes)
library(ggplot2)
library(dplyr)
library(janitor)

# Automatically use all available CPU cores (or up to your desired limit)
options(mc.cores = parallel::detectCores())

# Save compiled Stan models to disk so you don't recompile every run
rstan_options(auto_write = TRUE)

# Load Excel file into a DataFrame
#df <- read_excel("Y:/raw-imaging/Adrianna/experiments/analysis/Adrianna/NDNF/excels/static-patch_summary_data_rois_NDNF.xlsx", sheet = 1)
#df <- read_excel("Y:/raw-imaging/Adrianna/experiments/analysis/Adrianna/NDNF/excels/Natural-Images-4-repeats_summary_data_rois_NDNF.xlsx", sheet = 1)
df <- read_excel("Y:/raw-imaging/Adrianna/experiments/analysis/Adrianna/NDNF/excels/looming-stim_summary_data_rois_NDNF.xlsx", sheet = 1)

# Positive responsive rois
df_pos <- df[df$responsive == 1,]
df_pos$virus <- relevel(factor(df_pos$virus), ref = "sgRosa")
# Clean columns' names (remove special characters)
df_pos <- df_pos %>% 
  clean_names()

# Negative responsive rois
df_neg <- df[df$responsive == -1,]
df_neg$virus <- relevel(factor(df_neg$virus), ref = "sgRosa")
df_neg <- df_neg %>% 
  clean_names()

# print columns' names
names(df_pos)

priors <- c(
  prior(normal(0, 5), class = "b"),               # Effet du génotype
  prior(exponential(1), class = "sd"),            # Restreint la variance des individus/sessions
  prior(exponential(1), class = "sigma")          # Bruit résiduel des cellules
)

model <- brm(
  mean_avg_trial_win~ virus + (1 | mouse / session), 
  data = df_pos, #pos or neg responses
  family = gaussian(),
  chains = 4, cores = 4,
  iter = 4000,          # Augmenté pour résoudre le problème de Tail ESS
  warmup = 2000,
  prior = priors,
  control = list(
    adapt_delta = 0.999, # Force des pas d'échantillonnage plus petits (par défaut 0.80)
    max_treedepth = 15  # Évite un autre avertissement fréquent lié aux petits pas
  )
)

pairs(model, variable = c("sd_mouse__Intercept", "sd_mouse:session__Intercept", "sigma"))

summary(model)

# Affiche la distribution a posteriori et l'intervalle de crédibilité à 95%
posterior_summary(model, variable = "b_virussgCnr1") # Adaptez selon le nom exact de votre variable

# Alternativement, pour obtenir la probabilité que sgRosa > sgCnr1:
hyp <- hypothesis(model, "virussgCnr1 > 0")
hyp

post_prob <- hyp$hypothesis$Post.Prob
prob_pos  <- round(post_prob * 100)
prob_neg  <- 100 - prob_pos

label_pos <- paste0(prob_pos, "% Prob (> 0)")
label_neg <- paste0(prob_neg, "% Prob (< 0)")

# Extract MCMC draws for your treatment parameter
model %>%
  spread_draws(b_virussgCnr1) %>%
  ggplot(aes(x = b_virussgCnr1)) +
  stat_halfeye(
    aes(fill = after_stat(x > 0)),
    .width = c(0.90, 0.95),
    point_interval = mean_qi
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  scale_fill_manual(
    values = c("TRUE" = "#377eb8", "FALSE" = "#e41a1c"),
    labels = c("FALSE" = label_neg, "TRUE" = label_pos),
    name = "Posterior Probability"
  ) +
  labs(
    title = "Evoked Response Difference (1-Second Peak Window)",
    x = "Baseline-Subtracted ΔF/F Difference (KD - WT)",
    y = "Density"
  ) +
  theme_minimal()
