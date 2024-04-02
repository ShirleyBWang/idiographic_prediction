## U01 idiographic analyses 

##---- packages ----
library(tidyverse)
library(lubridate)
library(caret)
library(MLmetrics)
library(RANN)
library(psych)
library(reshape2)
library(corrplot)
library(tidymodels)
library(forecast)
library(conflicted)
library(emmeans)
library(lme4)
library(emmeans)
tidymodels_prefer()
#library(Metrics)

`%notin%` <- Negate(`%in%`)

##---- data ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# urge results
baseline_urge_ar1 <- read.csv("../results/baseline_ar1_urge_metrics.csv")
baseline_urge_affect <- read.csv("results/baseline_urge_affect.csv")
gp_urge_only <- read.csv("../results/gp_urge_only_window.csv")
gp_urge_affect <- read.csv("../results/urge_gp_affect_window.csv")

# merge urge results 
baseline_urge_ar1_tomerge <- baseline_urge_ar1 %>% 
  rename(mae_baseline_ar1 = mae, rmse_baseline_ar1 = rmse, r2_baseline_ar1 = r2)

baseline_urge_affect_tomerge <- baseline_urge_affect %>% 
  rename(mae_baseline_affect = mae, rmse_baseline_affect = rmse, r2_baseline_affect = r2)

gp_urge_only_tomerge <- gp_urge_only %>% 
  rename(mae_simple_gp = mae, rmse_simple_gp = rmse, r2_simple_gp = r2)

gp_urge_affect_tomerge <- gp_urge_affect %>% 
  rename(mae_gp_affect = mae, rmse_gp_affect = rmse, r2_gp_affect = r2)

urge_models_all <- baseline_urge_ar1_tomerge %>% 
  merge(baseline_urge_affect_tomerge, by = "ppt_id") %>% 
  merge(gp_urge_only_tomerge, by = "ppt_id") %>% 
  merge(gp_urge_affect_tomerge, by = "ppt_id")

# intent results
baseline_intent_ar1 <- read.csv("../results/baseline_ar1_intent_metrics.csv")
baseline_intent_affect <- read.csv("results/baseline_intent_affect.csv")
gp_intent_only <- read.csv("../results/gp_intent_only_window.csv")
gp_intent_affect <- read.csv("../results/intent_gp_affect_window.csv")

# merge intent results 
baseline_intent_ar1_tomerge <- baseline_intent_ar1 %>% 
  rename(mae_baseline_ar1 = mae, rmse_baseline_ar1 = rmse, r2_baseline_ar1 = r2)

baseline_intent_affect_tomerge <- baseline_intent_affect %>% 
  rename(mae_baseline_affect = mae, rmse_baseline_affect = rmse, r2_baseline_affect = r2)

gp_intent_only_tomerge <- gp_intent_only %>% 
  rename(mae_simple_gp = mae, rmse_simple_gp = rmse, r2_simple_gp = r2)

gp_intent_affect_tomerge <- gp_intent_affect %>% 
  rename(mae_gp_affect = mae, rmse_gp_affect = rmse, r2_gp_affect = r2)

intent_models_all <- baseline_intent_ar1_tomerge %>% 
  merge(baseline_intent_affect_tomerge, by = "ppt_id") %>% 
  merge(gp_intent_only_tomerge, by = "ppt_id") %>% 
  merge(gp_intent_affect_tomerge, by = "ppt_id")


##---- summary statistics ----
## urge
baseline_urge_ar1 %>% describe()
baseline_urge_affect %>% describe()
gp_urge_only %>% describe()
gp_urge_affect %>% describe()

## intent 
baseline_intent_ar1 %>% describe()
baseline_intent_affect %>% describe()
gp_intent_only %>% describe()
gp_intent_affect %>% describe()

##---- low, medium, large effect sizes ----
breaks <- c(0, 0.02, 0.13, 0.26, 1)

## urge
baseline_urge_ar1 %>% 
  mutate(interval = cut(r2, breaks = breaks, labels = c("<0.02", "0.02-0.13", "0.13-0.26", ">0.26"))) %>%
  group_by(interval) %>%
  summarize(count = n())

baseline_urge_affect %>% 
  mutate(interval = cut(r2, breaks = breaks, labels = c("<0.02", "0.02-0.13", "0.13-0.26", ">0.26"))) %>%
  group_by(interval) %>%
  summarize(count = n())

gp_urge_only %>% 
  mutate(interval = cut(r2, breaks = breaks, labels = c("<0.02", "0.02-0.13", "0.13-0.26", ">0.26"))) %>%
  group_by(interval) %>%
  summarize(count = n())

gp_urge_affect %>% 
  mutate(interval = cut(r2, breaks = breaks, labels = c("<0.02", "0.02-0.13", "0.13-0.26", ">0.26"))) %>%
  group_by(interval) %>%
  summarize(count = n())

## intent
baseline_intent_ar1 %>% 
  mutate(interval = cut(r2, breaks = breaks, labels = c("<0.02", "0.02-0.13", "0.13-0.26", ">0.26"))) %>%
  group_by(interval) %>%
  summarize(count = n())

baseline_intent_affect %>% 
  mutate(interval = cut(r2, breaks = breaks, labels = c("<0.02", "0.02-0.13", "0.13-0.26", ">0.26"))) %>%
  group_by(interval) %>%
  summarize(count = n())

gp_intent_only %>% 
  mutate(interval = cut(r2, breaks = breaks, labels = c("<0.02", "0.02-0.13", "0.13-0.26", ">0.26"))) %>%
  group_by(interval) %>%
  summarize(count = n())

gp_intent_affect %>% 
  mutate(interval = cut(r2, breaks = breaks, labels = c("<0.02", "0.02-0.13", "0.13-0.26", ">0.26"))) %>%
  group_by(interval) %>%
  summarize(count = n())

##---- statistical comparison ----
urge_r2_melted <- urge_models_all %>% 
  select(contains("r2"), ppt_id) %>% melt(id.vars = "ppt_id") 
urge_r2_lmer_null <- lmerTest::lmer(value ~ (1 | ppt_id), data = urge_r2_melted)
urge_r2_lmer <- lmerTest::lmer(value ~ variable + (1 | ppt_id), data = urge_r2_melted)
anova(urge_r2_lmer_null, urge_r2_lmer)
lsmeans::lsmeans(urge_r2_lmer, pairwise ~ variable, adjust = "tukey")

intent_r2_melted <- intent_models_all %>% 
  select(contains("r2"), ppt_id) %>% melt(id.vars = "ppt_id") 
intent_r2_lmer_null <- lmerTest::lmer(value ~ (1 | ppt_id), data = intent_r2_melted)
intent_r2_lmer <- lmerTest::lmer(value ~ variable + (1 | ppt_id), data = intent_r2_melted)
anova(intent_r2_lmer_null, intent_r2_lmer)
lsmeans::lsmeans(intent_r2_lmer, pairwise ~ variable, adjust = "tukey")

##---- prep for plot ----
## set RMSE of the high intent (affect gp) to NA
which(gp_intent_affect$rmse > 30)
gp_intent_affect$rmse[67] <- NA

## set RMSE of the high urge (affect gp) to NA
which(gp_urge_affect$rmse > 10)
gp_urge_affect$rmse[39] <- NA

##---- violin plots for intent ----
baseline_ar1_intent_melted <- melt(baseline_intent_ar1, value.name = "Baseline (AR1)", variable.name = "Metric")
baseline_intent_affect_melted <- melt(baseline_intent_affect, value.name = "BaselineAffect", variable.name = "Metric")
gp_intent_melted <- melt(gp_intent_only, value.name = "SimpleGP", variable.name = "Metric")
gp_intent_affect_melted <- melt(gp_intent_affect, value.name = "GPAffect", variable.name = "Metric")

#baseline_intent_ar1$ppt_id[which(baseline_intent_ar1$ppt_id %notin% baseline_intent_affect$ppt_id)]

intent_violin_dat <- as.data.frame(cbind(baseline_ar1_intent_melted, 
                                         baseline_intent_affect_melted$BaselineAffect,
                                         gp_intent_melted$SimpleGP,
                                         gp_intent_affect_melted$GPAffect))
names(intent_violin_dat) <- c("ppt_id", "Metric", "Baseline \n(AR1)", "Baseline \n(Affect)", "GP \n(Simple)", "GP \n(Affect)")
intent_violin_dat$Metric <- recode(intent_violin_dat$Metric, "mae"="MAE", "rmse"="RMSE", "r2" = "R2")
intent_violin_melted <- melt(intent_violin_dat)

intent_violin_melted$Metric <- intent_violin_melted$Metric %>% 
  fct_relevel(., "R2", "RMSE", "MAE")

## ggplot
intent_violin_melted %>% 
  ggplot(., aes(x = variable, y = value, color = after_scale(alpha(fill, 0.6)), fill = variable)) +
  geom_point(size = 0.7) +
  geom_violin(size = 0.6, alpha = 0.1) +
  geom_boxplot(width=0.1, alpha = 0.2, fatten = 3, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 6) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.15) +
  ggtitle("Idiographic Prediction Model Metrics: \nPredicting Suicidal Intent") +
  facet_wrap(~Metric, scales = "free") +
  theme_bw() + 
  theme(axis.text=element_text(size = 10), axis.title=element_blank(), 
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = "none", 
        strip.text = element_text(size = 10)) +
  guides(fill = FALSE) +
  scale_color_manual(values = c("#8B8B8B", "#789B60", "#56A2A5", "#6F4988")) +
  scale_fill_manual(values = c("#8B8B8B", "#789B60", "#56A2A5", "#6F4988"))
ggsave("figs/intent_mod_comparison.pdf", bg='transparent', width = 12, height = 6)

##---- violin plots for urge ----
baseline_ar1_urge_melted <- melt(baseline_urge_ar1, value.name = "Baseline (AR1)", variable.name = "Metric")
baseline_urge_affect_melted <- melt(baseline_urge_affect, value.name = "BaselineAffect", variable.name = "Metric")
gp_urge_melted <- melt(gp_urge_only, value.name = "SimpleGP", variable.name = "Metric")
gp_urge_affect_melted <- melt(gp_urge_affect, value.name = "GPAffect", variable.name = "Metric")

urge_violin_dat <- as.data.frame(cbind(baseline_ar1_urge_melted, 
                                       baseline_urge_affect_melted$BaselineAffect,
                                       gp_urge_melted$SimpleGP,
                                       gp_urge_affect_melted$GPAffect))
names(urge_violin_dat) <- c("ppt_id", "Metric", "Baseline \n(AR1)", "Baseline \n(Affect)", "GP \n(Simple)", "GP \n(Affect)")
urge_violin_dat$Metric <- recode(urge_violin_dat$Metric, "mae"="MAE", "rmse"="RMSE", "r2" = "R2")
urge_violin_melted <- melt(urge_violin_dat)

urge_violin_melted$Metric <- urge_violin_melted$Metric %>% 
  fct_relevel(., "R2", "RMSE", "MAE")

## ggplot
urge_violin_melted %>% 
  ggplot(., aes(x = variable, y = value, color = after_scale(alpha(fill, 0.6)), fill = variable)) +
  geom_point(size = 0.7) +
  geom_violin(size = 0.6, alpha = 0.1) +
  geom_boxplot(width=0.1, alpha = 0.2, fatten = 3, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 6) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.15) +
  ggtitle("Idiographic Prediction Model Metrics: \nPredicting Suicidal Urges") +
  facet_wrap(~Metric, scales = "free") +
  theme_bw() + 
  theme(axis.text=element_text(size = 10), axis.title=element_blank(), 
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = "none", 
        strip.text = element_text(size = 10)) +
  guides(fill = FALSE) +
  scale_color_manual(values = c("#8B8B8B", "#789B60", "#56A2A5", "#6F4988")) +
  scale_fill_manual(values = c("#8B8B8B", "#789B60", "#56A2A5", "#6F4988"))
ggsave("figs/urge_mod_comparison.pdf", bg='transparent', width = 12, height = 6)


##---- Urge correlations ----
u01_urge <- read.csv("data/urge_affect_lagged_021623.csv")
u01_urge$start_time <- ymd_hms(u01_urge$start_time)

# calculate features
urge_ft <- extract_features(u01_urge, group_var = "ppt_id", value_var = "sitb_si_urge", features = "all")
urge_ft_df <- features_to_df(urge_ft, data.format="wide", group_var = "ppt_id") %>% 
  select(ppt_id, f.min_sitb_si_urge, f.mean_sitb_si_urge, f.max_sitb_si_urge, f.sd_sitb_si_urge, 
         f.mssd.lag1_sitb_si_urge,f.perc.unique.values_sitb_si_urge,
         f.max.change_sitb_si_urge, f.prob.acute.change.9_sitb_si_urge)

urge_comp <- merge(urge_models_all, urge_ft_df, by = "ppt_id") %>% 
  select(-starts_with(c("mae", "rmse"))) %>% 
  select(-ppt_id) %>% 
  rename("R2 Baseline (ar1)" = r2_baseline_ar1, 
         "R2 Baseline (affect)" = r2_baseline_affect,
         "R2 GP (simple)" = r2_simple_gp,
         "R2 GP (affect)" = r2_gp_affect, 
         "Urge Minimum" = f.min_sitb_si_urge, 
         "Urge Mean" = f.mean_sitb_si_urge, 
         "Urge Maximum" = f.max_sitb_si_urge, 
         "Urge SD" = f.sd_sitb_si_urge, 
         "Urge MSSD" = f.mssd.lag1_sitb_si_urge, 
         "Urge Perc Unique Vals" = f.perc.unique.values_sitb_si_urge,
         "Urge Maximum Change" = f.max.change_sitb_si_urge, 
         "Urge Prob Acute Change" = f.prob.acute.change.9_sitb_si_urge)

urge_cormat <- urge_comp %>% 
  cor(., use = "pairwise.complete.obs")

res1 <- cor.mtest(urge_comp, conf.level = .95)

pdf("figs/urge_models/urge_correlations.pdf", width = 10, height = 10)
corrplot(urge_cormat, method = "color", 
         #addCoef.col = "gray40", 
         insig = "blank", 
         tl.col = "black", 
         type = "lower",
         number.cex = 1.2, 
        # p.mat = res1$p, 
         bg = "lightgrey",
         pch.cex = 1.2, cl.pos = "b")
dev.off()

##---- Intent correlations ----
u01_intent <- read.csv("data/intent_affect_lagged_021623.csv")
u01_intent$start_time <- ymd_hms(u01_intent$start_time)

# calculate features
intent_ft <- extract_features(u01_intent, group_var = "ppt_id", value_var = "sitb_si_intent", features = "all")
intent_ft_df <- features_to_df(intent_ft, data.format="wide", group_var = "ppt_id") %>% 
  select(ppt_id, f.min_sitb_si_intent, f.mean_sitb_si_intent, f.max_sitb_si_intent, f.sd_sitb_si_intent, 
         f.mssd.lag1_sitb_si_intent,f.perc.unique.values_sitb_si_intent,
         f.max.change_sitb_si_intent, f.prob.acute.change.9_sitb_si_intent)

intent_comp <- merge(intent_models_all, intent_ft_df, by = "ppt_id") %>% 
  select(-starts_with(c("mae", "rmse"))) %>% 
  select(-ppt_id) %>% 
  rename("R2 Baseline (ar1)" = r2_baseline_ar1, 
         "R2 Baseline (affect)" = r2_baseline_affect,
         "R2 GP (simple)" = r2_simple_gp,
         "R2 GP (affect)" = r2_gp_affect, 
         "Intent Minimum" = f.min_sitb_si_intent, 
         "Intent Mean" = f.mean_sitb_si_intent, 
         "Intent Maximum" = f.max_sitb_si_intent, 
         "Intent SD" = f.sd_sitb_si_intent, 
         "Intent MSSD" = f.mssd.lag1_sitb_si_intent, 
         "Intent Perc Unique Vals" = f.perc.unique.values_sitb_si_intent,
         "Intent Maximum Change" = f.max.change_sitb_si_intent, 
         "Intent Prob Acute Change" = f.prob.acute.change.9_sitb_si_intent)

intent_cormat <- intent_comp %>% 
  cor(., use = "pairwise.complete.obs")

res1 <- cor.mtest(intent_comp, conf.level = .95)

pdf("figs/intent_models/intent_correlations.pdf", width = 10, height = 10)
corrplot(intent_cormat, method = "color", 
         #addCoef.col = "gray40", 
         insig = "blank", 
         tl.col = "black", 
         type = "lower",
         number.cex = 1.2, 
         #p.mat = res1$p, 
         bg = "lightgrey",
         pch.cex = 1.2, cl.pos = "b")
dev.off()





