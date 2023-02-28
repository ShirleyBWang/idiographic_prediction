## U01 idiographic analyses - for diss

##---- packages ----
library(tidyverse)
library(lubridate)
library(caret)
library(MLmetrics)
library(RANN)
library(tsfeaturex)
library(glmnet)
library(hydroGOF) 
library(psych)
library(reshape2)
library(corrplot)
library(tidymodels)
library(forecast)
library(conflicted)
tidymodels_prefer()
conflict_prefer("rmse", "yardstick")
#library(Metrics)

`%notin%` <- Negate(`%in%`)

## ggplot theme for papers/talks
paper_color = 'black'
paper_theme <- theme(
  axis.text = element_text(size = 11, color = paper_color), 
  text = element_text(color = paper_color),
  plot.title = element_text(size = 15, hjust = 0.5),
  strip.text = element_text(size = 11, color = paper_color),
  panel.background = element_rect(fill='transparent'), #transparent panel bg
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  panel.grid.major = element_blank(), #remove major gridlines
  panel.grid.minor = element_blank(), #remove minor gridlines
  legend.background = element_rect(fill='transparent'), #transparent legend bg
  legend.box.background = element_rect(fill='transparent') #transparent legend panel
)

##---- data ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# urge
u01_urge <- read.csv("data/urge_affect_lagged_021623.csv")
u01_urge$start_time <- ymd_hms(u01_urge$start_time)

affect_vars <- u01_urge %>% 
  select(starts_with("affect_")) %>% 
  select(-c(affect_desire_approach, affect_desire_avoid)) %>%
  select(-starts_with("affect_stress_")) %>%
  names()

u01_urge %>% select(affect_vars) %>% describe() 
# remove impulsive, overwhelmed (total n ~ 45,000 compared to 62,000 full responses)
# consider removing happy, relaxed, sad, stressed, tense 

affect_vars <- u01_urge %>% 
  select(starts_with("affect_")) %>% 
  select(-c(affect_desire_approach, affect_desire_avoid)) %>%
  select(-starts_with("affect_stress_")) %>%
  select(-ends_with(c("impulsive", "overwhelmed", "happy", "relaxed", "sad", "stressed", "tense"))) %>% 
  names()

u01_urge %>% select(affect_vars) %>% describe() 
u01_urge %>% group_by(ppt_id) %>% count() # 88 people 

## intent
u01_intent <- read.csv("data/intent_affect_lagged_021623.csv")
u01_intent$start_time <- ymd_hms(u01_intent$start_time)
u01_intent %>% group_by(ppt_id) %>% count() # 89 people 

##---- descriptives ----
unique(u01_intent$ppt_id) %>% 
  str_detect("mgh") %>% 
  summary() # 43 adults, 42 adol

# how many surveys did people complete? 
u01_intent %>% 
  group_by(ppt_id) %>% 
  count() %>% 
  describe()

# how long were people responding for? 
u01_intent %>% 
  group_by(ppt_id) %>% 
  summarise(days_in_study = max(day_in_study)) %>% 
  describe()

# what was compliance? 
compliance <- u01_intent %>% 
  group_by(ppt_id) %>% 
  summarise(tot_surveys = max(day_in_study) * 6,
            complete_surveys = n()) %>% 
  mutate(compliance = complete_surveys / tot_surveys)

compliance <- 
  compliance %>% 
  mutate(compliance_adjusted = ifelse(compliance > 1, 1, compliance))
compliance$compliance_adjusted %>% describe()

# how long between prompts? 
u01_intent$time_lead %>% describe()
u01_urge$time_lead %>% describe()

##---- Functions: AR1 Urge Models----
ar_model_baseline_urge <- function(df, train_perc) {
  
  # initiate vectors to store information 
  n_obs <- nrow(df)
  time_in_study <- vector()
  y_true <- vector()
  y_pred <- vector()
  
  # get training and test data 
  train_num <- round(train_perc * nrow(df)) # get nrows for train
  test_num <- nrow(df) - train_num # get nrows for test
  
  train_dat <- df %>% slice_head(n = train_num) %>% 
    select(ppt_id, urge_lead, sitb_si_urge, time_lead, time_in_study)
  test_dat <- df %>% slice_tail(n = test_num) %>% 
    select(ppt_id, urge_lead, sitb_si_urge, time_lead, time_in_study)
  
  # save initial values
  time_in_study <- df$time_in_study
  y_true <- df$urge_lead
  y_pred[1:train_num] <- NA
  
  # set model specifications
  lm_mod <- linear_reg() %>% set_engine("lm") %>% set_mode("regression")
  
  # build model, make predictions
  ar1_fit <- lm_mod %>% 
    fit(urge_lead ~ sitb_si_urge, data = train_dat)
  
  ar1_pred <- predict(ar1_fit, new_data = test_dat)
  
  # save results
  y_pred <- append(y_pred, ar1_pred$.pred)

  ar_results <- tibble(time_in_study, y_true, y_pred)
  mae_val <- yardstick::mae(ar_results, y_true, y_pred)$.estimate
  rmse_val <- yardstick::rmse(ar_results, y_true, y_pred)$.estimate
  r2_val <- yardstick::rsq(ar_results, y_true, y_pred)$.estimate
  ar_metrics <- tibble(df$ppt_id[1], mae_val, rmse_val, r2_val)
  names(ar_metrics) <- c("ppt_id", "mae", "rmse", "r2")
  return(list("preds" = ar_results, "metrics" = ar_metrics))
}

plot_ar1_urge <- function(ar1_output, ppt_df, display = FALSE) {
  ## restructure data for plotting (lag y_pred and y_true back by one)
  preds <- ar1_output$preds %>% 
    mutate(pred_urge = dplyr::lag(y_pred),
           sitb_si_urge = ppt_df$sitb_si_urge)
  
  ## create ggplot
  ar1plot <- ggplot(preds, aes(x = time_in_study, y = sitb_si_urge)) +
    geom_line(col = '#CC6633') +
    geom_point(col = '#CC6633', size = 0.75) +
    geom_line(aes(y = pred_urge), col = '#006666') + 
    geom_point(aes(y = pred_urge), col = '#006666', size = 0.75) +
    labs(x = 'hours of study', y = 'Suicidal Urges Severity', 
         title = sprintf("Baseline Model (AR-1) Predicting Suicidal Urges for %s \nmae = %g, rmse = %g, r2 = %g", 
                         ar1_output$metrics$ppt_id, round(ar1_output$metrics$mae, 3), 
                         round(ar1_output$metrics$rmse, 3), round(ar1_output$metrics$r2, 3)), font.main = 1) +
    ylim(0, 10) + 
    scale_x_continuous() +
    paper_theme
  
  if (display == TRUE){
    show(ar1plot)
  }
  
  ## save ggplot
  ggsave(sprintf("figs/urge_models/baseline_urge_ar1/urge_baseline_ar1_%s.pdf", ar1_output$metrics$ppt_id), 
         bg='transparent', width = 8, height = 4)
}

## try on just one patient (fc003)
fc025 <- u01_urge %>% filter(ppt_id == 'u01-fc025')
ar_fc025 <- ar_model_baseline_urge(df = fc025, train_perc = 0.75)
ar_fc025$metrics
plot_ar1_urge(ar_fc025, fc025, display = TRUE)

##---- Fit AR1 Urge Models----
ar1_metrics_urge <- tibble(ppt_id = character(),
                      mae = numeric(),
                      rmse = numeric(), 
                      r2 = numeric())
#ar1_metrics <- ar1_metrics %>% add_row(ar_fc003$metrics)
for (ppt in unique(u01_urge$ppt_id)) {
  # subset data
  this_ppt <- u01_urge %>% filter(ppt_id == ppt)
  this_mod <- ar_model_baseline_urge(df = this_ppt, train_perc = 0.5)
  ar1_metrics_urge <- ar1_metrics_urge %>% add_row(this_mod$metrics)
  plot_ar1_urge(this_mod, this_ppt, display = FALSE)
  
  # keep track of progress
  print(sprintf("%s: mae = %g, rmse = %g, r2 = %g", 
                this_mod$metrics$ppt_id, round(this_mod$metrics$mae, 3), 
                round(this_mod$metrics$rmse, 3), round(this_mod$metrics$r2, 3)))
  
}

summary(ar1_metrics_urge)
write.csv(ar1_metrics_urge, "results/baseline_ar1_urge_metrics.csv", row.names = FALSE)

##---- Functions: AR1 Intent Models ----
ar_model_baseline_intent <- function(df, train_perc) {
  
  # initiate vectors to store information 
  n_obs <- nrow(df)
  time_in_study <- vector()
  y_true <- vector()
  y_pred <- vector()
  
  # get training and test data 
  train_num <- round(train_perc * nrow(df)) # get nrows for train
  test_num <- nrow(df) - train_num # get nrows for test
  
  train_dat <- df %>% slice_head(n = train_num) %>% 
    select(ppt_id, intent_lead, sitb_si_intent, time_lead, time_in_study)
  test_dat <- df %>% slice_tail(n = test_num) %>% 
    select(ppt_id, intent_lead, sitb_si_intent, time_lead, time_in_study)
  
  # save initial values
  time_in_study <- df$time_in_study
  y_true <- df$intent_lead
  y_pred[1:train_num] <- NA
  
  # set model specifications
  lm_mod <- linear_reg() %>% set_engine("lm") %>% set_mode("regression")
  
  # build model, make predictions
  ar1_fit <- lm_mod %>% 
    fit(intent_lead ~ sitb_si_intent, data = train_dat)
  
  ar1_pred <- predict(ar1_fit, new_data = test_dat)
  
  # save results
  y_pred <- append(y_pred, ar1_pred$.pred)
  ar_results <- tibble(time_in_study, y_true, y_pred)
  mae_val <- yardstick::mae(ar_results, y_true, y_pred)$.estimate
  rmse_val <- yardstick::rmse(ar_results, y_true, y_pred)$.estimate
  r2_val <- yardstick::rsq(ar_results, y_true, y_pred)$.estimate
  ar_metrics <- tibble(df$ppt_id[1], mae_val, rmse_val, r2_val)
  names(ar_metrics) <- c("ppt_id", "mae", "rmse", "r2")
  return(list("preds" = ar_results, "metrics" = ar_metrics))
}

plot_ar1_intent <- function(ar1_output, ppt_df, display = FALSE) {
  ## restructure data for plotting (lag y_pred and y_true back by one)
  preds <- ar1_output$preds %>% 
    mutate(pred_intent = dplyr::lag(y_pred),
           sitb_si_intent = ppt_df$sitb_si_intent)
  
  ## create ggplot
  ar1plot <- ggplot(preds, aes(x = time_in_study, y = sitb_si_intent)) +
    geom_line(col = '#CC6633') +
    geom_point(col = '#CC6633', size = 0.75) +
    geom_line(aes(y = pred_intent), col = '#006666') + 
    geom_point(aes(y = pred_intent), col = '#006666', size = 0.75) +
    labs(x = 'hours of study', y = 'Suicidal Intent Severity', 
         title = sprintf("Baseline Model (AR-1) Predicting Suicidal Intent for %s \nmae = %g, rmse = %g, r2 = %g", 
                         ar1_output$metrics$ppt_id, round(ar1_output$metrics$mae, 3), 
                         round(ar1_output$metrics$rmse, 3), round(ar1_output$metrics$r2, 3)), font.main = 1) +
    ylim(0, 10) + 
    scale_x_continuous() +
    paper_theme
  
  if (display == TRUE){
    show(ar1plot)
  }
  
  ## save ggplot
  ggsave(sprintf("figs/intent_models/baseline_intent_ar1/intent_baseline_ar1_%s.pdf", ar1_output$metrics$ppt_id), 
         bg='transparent', width = 8, height = 4)
}

## try on just one patient (fc003)
fc010 <- u01_intent %>% filter(ppt_id == 'u01-fc010')
ar_fc010 <- ar_model_baseline_intent(df = fc010, train_perc = 0.75)
ar_fc010$metrics
plot_ar1_intent(ar_fc010, fc010, display = TRUE)


##---- Fit AR1 Intent Models ----
ar1_metrics_intent <- tibble(ppt_id = character(),
                             mae = numeric(),
                             rmse = numeric(), 
                             r2 = numeric())

for (ppt in unique(u01_intent$ppt_id)) {
  # subset data
  this_ppt <- u01_intent %>% filter(ppt_id == ppt)
  this_mod <- ar_model_baseline_intent(df = this_ppt, train_perc = 0.5)
  ar1_metrics_intent <- ar1_metrics_intent %>% add_row(this_mod$metrics)
  plot_ar1_intent(this_mod, this_ppt, display = FALSE)
  
  # keep track of progress
  print(sprintf("%s: mae = %g, rmse = %g, r2 = %g", 
                this_mod$metrics$ppt_id, round(this_mod$metrics$mae, 3), 
                round(this_mod$metrics$rmse, 3), round(this_mod$metrics$r2, 3)))
  
}

summary(ar1_metrics_intent)
write.csv(ar1_metrics_intent, "results/baseline_ar1_intent_metrics.csv", row.names = FALSE)

##---- Functions: Urge Model with Affect ----
affect_model_baseline_urge <- function(df, train_perc) {
  
  this_ppt_id = df$ppt_id[1]
  
  # remove missing rows 
  #df <- df %>% 
  #  select(urge_lead, sitb_si_urge, affect_vars) %>% 
  #  na.omit()
  
  # initiate vectors to store information 
  n_obs <- nrow(df)
  time_in_study <- vector()
  y_true <- vector()
  y_pred <- vector()
  
  # get training and test data 
  train_num <- round(train_perc * nrow(df)) # get nrows for train
  test_num <- nrow(df) - train_num # get nrows for test
  
  train_dat <- df %>% slice_head(n = train_num) %>% 
    select(urge_lead, sitb_si_urge, affect_vars)
  test_dat <- df %>% slice_tail(n = test_num) %>% 
    select(urge_lead, sitb_si_urge, affect_vars)
  
  if (is.constant(na.omit(train_dat$urge_lead))) {
    print("skipping this patient, constant urge_lead")
    return(0)
  } else if (sum(train_dat$urge_lead != 0) < 2) {
    print("skipping this patient, near-constant urge_lead")
    return(0)
  } else {
    # save initial values
    time_in_study <- df$time_in_study
    y_true <- df$urge_lead
    y_pred[1:train_num] <- NA
    
    # set up 5-fold CV
    set.seed(123)
    urge_folds <- vfold_cv(data = train_dat, v = 5, repeats = 3, strata = urge_lead)
    
    # set model specifications
    glmnet_model <- linear_reg(penalty = tune(), mixture = tune()) %>% 
      set_engine("glmnet") %>% 
      set_mode("regression")
    
    # set recipe
    urge_recipe <- 
      recipe(urge_lead ~ ., data = train_dat) %>% 
      step_impute_knn(all_predictors()) %>%
      step_nzv(all_predictors()) %>% 
      step_corr(all_predictors()) %>%
      step_lincomb(all_predictors()) %>% 
      step_normalize(all_predictors())
      
    
    # workflow and metrics
    urge_wflow <- 
      workflow() %>% 
      add_model(glmnet_model) %>% 
      add_recipe(urge_recipe)
    
    urge_ms <- metric_set(mae, rmse, rsq)
    
    # tune model and select best hyperparameters
    urge_tune <- 
      urge_wflow %>%
      tune_grid(
        resamples = urge_folds, 
        grid = 5,
        metrics = urge_ms
      )
    urge_param_final <- select_best(urge_tune, metric = "rmse")
    
    # finalize workflow with the best hyperparameter values
    urge_wflow_final <- 
      urge_wflow %>%
      finalize_workflow(urge_param_final)
  
    # fit and predict
    urge_deploy <- 
      urge_wflow_final %>%
      fit(train_dat)
    
    urge_predict <- predict(urge_deploy, new_data = test_dat)
    
    # save results
    y_pred <- append(y_pred, urge_predict$.pred)
    
    ar_results <- tibble(time_in_study, y_true, y_pred)
    mae_val <- yardstick::mae(ar_results, y_true, y_pred)$.estimate
    rmse_val <- yardstick::rmse(ar_results, y_true, y_pred)$.estimate
    r2_val <- yardstick::rsq(ar_results, y_true, y_pred)$.estimate
    ar_metrics <- tibble(this_ppt_id, mae_val, rmse_val, r2_val)
    names(ar_metrics) <- c("ppt_id", "mae", "rmse", "r2")
    return(list("preds" = ar_results, "metrics" = ar_metrics))
  }
}

plot_affect_urge <- function(affect_output, ppt_df, display = FALSE) {
  ## restructure data for plotting (lag y_pred and y_true back by one)
  preds <- affect_output$preds %>% 
    mutate(pred_urge = dplyr::lag(y_pred),
           sitb_si_urge = ppt_df$sitb_si_urge)
  
  ## create ggplot
  ar1plot <- ggplot(preds, aes(x = time_in_study, y = sitb_si_urge)) +
    geom_line(col = '#CC6633') +
    geom_point(col = '#CC6633', size = 0.75) +
    geom_line(aes(y = pred_urge), col = '#006666') + 
    geom_point(aes(y = pred_urge), col = '#006666', size = 0.75) +
    labs(x = 'hours of study', y = 'Suicidal Urges Severity', 
         title = sprintf("Baseline Model (with Affect) Predicting Suicidal Urges for %s \nmae = %g, rmse = %g, r2 = %g", 
                         affect_output$metrics$ppt_id, round(affect_output$metrics$mae, 3), 
                         round(affect_output$metrics$rmse, 3), round(affect_output$metrics$r2, 3)), font.main = 1) +
    ylim(0, 10) + 
    scale_x_continuous() +
    paper_theme
  
  if (display == TRUE){
    show(ar1plot)
  }
  
  ## save ggplot
  ggsave(sprintf("figs/urge_models/baseline_urge_affect/urge_baseline_affeect_%s.pdf", affect_output$metrics$ppt_id), 
         bg='transparent', width = 8, height = 4)
}

## try on just one patient (mgh264)
mgh267 <- u01_urge %>% filter(ppt_id == 'u01-mgh267')
affect_mgh267 <- affect_model_baseline_urge(df = mgh267, train_perc = 0.5)
affect_mgh267$metrics %>%
plot_affect_urge(affect_mgh267, mgh267, display = FALSE)

##---- Fit Urge Models with Affect ----
affect_metrics_urge <- tibble(ppt_id = character(),
                              mae = numeric(),
                              rmse = numeric(), 
                              r2 = numeric())
#ar1_metrics <- ar1_metrics %>% add_row(ar_fc003$metrics)
for (ppt in unique(u01_urge$ppt_id)) {
  # subset data
  this_ppt <- u01_urge %>% filter(ppt_id == ppt)
  this_mod <- affect_model_baseline_urge(df = this_ppt, train_perc = 0.5)
  
  if(is.list(this_mod)){ 
  
    affect_metrics_urge <- affect_metrics_urge %>% add_row(this_mod$metrics)
    plot_affect_urge(this_mod, this_ppt, display = FALSE)
    # keep track of progress
    print(sprintf("%s: mae = %g, rmse = %g, r2 = %g", 
                  this_mod$metrics$ppt_id, round(this_mod$metrics$mae, 3), 
                  round(this_mod$metrics$rmse, 3), round(this_mod$metrics$r2, 3)))
  }
  else {
    NA_metrics <- tibble(ppt_id = ppt,
                         mae = NA,
                         rmse = NA,
                         r2 = NA)
    affect_metrics_urge <- affect_metrics_urge %>% add_row(NA_metrics)
  }
  
}

summary(affect_metrics_urge)
write.csv(affect_metrics_urge, "results/baseline_urge_affect.csv", row.names = FALSE)
summary(ar1_metrics_urge) # compare to AR1


##---- Functions: Intent Model with Affect ----
affect_model_baseline_intent <- function(df, train_perc) {
  
  this_ppt_id = df$ppt_id[1]
  
  # remove missing rows 
  #df <- df %>% 
  #  select(intent_lead, sitb_si_intent, affect_vars) %>% 
  #  na.omit()
  
  # initiate vectors to store information 
  n_obs <- nrow(df)
  time_in_study <- vector()
  y_true <- vector()
  y_pred <- vector()
  
  # get training and test data 
  train_num <- round(train_perc * nrow(df)) # get nrows for train
  test_num <- nrow(df) - train_num # get nrows for test
  
  train_dat <- df %>% slice_head(n = train_num) %>% 
    select(intent_lead, sitb_si_intent, affect_vars)
  test_dat <- df %>% slice_tail(n = test_num) %>% 
    select(intent_lead, sitb_si_intent, affect_vars)
  
  if (is.constant(na.omit(train_dat$intent_lead))) {
    print("skipping this patient, constant intent_lead")
    return(0)
  } else if (sum(train_dat$intent_lead != 0) < 2) {
    print("skipping this patient, near-constant intent_lead")
    return(0)
  } else {
    # save initial values
    time_in_study <- df$time_in_study
    y_true <- df$intent_lead
    y_pred[1:train_num] <- NA
    
    # set up 5-fold CV
    set.seed(123)
    intent_folds <- vfold_cv(data = train_dat, v = 5, repeats = 3, strata = intent_lead)
    
    # set model specifications
    glmnet_model <- linear_reg(penalty = tune(), mixture = tune()) %>% 
      set_engine("glmnet") %>% 
      set_mode("regression")
    
    # set recipe
    intent_recipe <- 
      recipe(intent_lead ~ ., data = train_dat) %>% 
      step_impute_knn(all_predictors()) %>%
      step_nzv(all_predictors()) %>% 
      step_corr(all_predictors()) %>%
      step_lincomb(all_predictors()) %>% 
      step_normalize(all_predictors())
    
    
    # workflow and metrics
    intent_wflow <- 
      workflow() %>% 
      add_model(glmnet_model) %>% 
      add_recipe(intent_recipe)
    
    intent_ms <- metric_set(mae, rmse, rsq)
    
    # tune model and select best hyperparameters
    intent_tune <- 
      intent_wflow %>%
      tune_grid(
        resamples = intent_folds, 
        grid = 5,
        metrics = intent_ms
      )
    intent_param_final <- select_best(intent_tune, metric = "rmse")
    
    # finalize workflow with the best hyperparameter values
    intent_wflow_final <- 
      intent_wflow %>%
      finalize_workflow(intent_param_final)
    
    # fit and predict
    intent_deploy <- 
      intent_wflow_final %>%
      fit(train_dat)
    
    intent_predict <- predict(intent_deploy, new_data = test_dat)
    
    # save results
    y_pred <- append(y_pred, intent_predict$.pred)
    
    ar_results <- tibble(time_in_study, y_true, y_pred)
    mae_val <- yardstick::mae(ar_results, y_true, y_pred)$.estimate
    rmse_val <- yardstick::rmse(ar_results, y_true, y_pred)$.estimate
    r2_val <- yardstick::rsq(ar_results, y_true, y_pred)$.estimate
    ar_metrics <- tibble(this_ppt_id, mae_val, rmse_val, r2_val)
    names(ar_metrics) <- c("ppt_id", "mae", "rmse", "r2")
    return(list("preds" = ar_results, "metrics" = ar_metrics))
  }
}

plot_affect_intent <- function(affect_output, ppt_df, display = FALSE) {
  ## restructure data for plotting (lag y_pred and y_true back by one)
  preds <- affect_output$preds %>% 
    mutate(pred_intent = dplyr::lag(y_pred),
           sitb_si_intent = ppt_df$sitb_si_intent)
  
  ## create ggplot
  ar1plot <- ggplot(preds, aes(x = time_in_study, y = sitb_si_intent)) +
    geom_line(col = '#CC6633') +
    geom_point(col = '#CC6633', size = 0.75) +
    geom_line(aes(y = pred_intent), col = '#006666') + 
    geom_point(aes(y = pred_intent), col = '#006666', size = 0.75) +
    labs(x = 'hours of study', y = 'Suicidal Intent Severity', 
         title = sprintf("Baseline Model (with Affect) Predicting Suicidal Intent for %s \nmae = %g, rmse = %g, r2 = %g", 
                         affect_output$metrics$ppt_id, round(affect_output$metrics$mae, 3), 
                         round(affect_output$metrics$rmse, 3), round(affect_output$metrics$r2, 3)), font.main = 1) +
    ylim(0, 10) + 
    scale_x_continuous() +
    paper_theme
  
  if (display == TRUE){
    show(ar1plot)
  }
  
  ## save ggplot
  ggsave(sprintf("figs/intent_models/baseline_intent_affect/intent_baseline_affeect_%s.pdf", affect_output$metrics$ppt_id), 
         bg='transparent', width = 8, height = 4)
}

## try on just one patient (mgh264)
mgh268 <- u01_intent %>% filter(ppt_id == 'u01-mgh268')
affect_mgh268 <- affect_model_baseline_intent(df = mgh268, train_perc = 0.5)
affect_mgh268$metrics %>% View()
plot_affect_intent(affect_mgh268, mgh268, display = FALSE)

##---- Fit Intent Models with Affect ----
affect_metrics_intent <- tibble(ppt_id = character(),
                                mae = numeric(),
                                rmse = numeric(), 
                                r2 = numeric())
#ar1_metrics <- ar1_metrics %>% add_row(ar_fc003$metrics)
for (ppt in unique(u01_intent$ppt_id)) {
  # subset data
  this_ppt <- u01_intent %>% filter(ppt_id == ppt)
  this_mod <- affect_model_baseline_intent(df = this_ppt, train_perc = 0.5)
  
  if(is.list(this_mod)){ 
    
    affect_metrics_intent <- affect_metrics_intent %>% add_row(this_mod$metrics)
    plot_affect_intent(this_mod, this_ppt, display = FALSE)
    # keep track of progress
    print(sprintf("%s: mae = %g, rmse = %g, r2 = %g", 
                  this_mod$metrics$ppt_id, round(this_mod$metrics$mae, 3), 
                  round(this_mod$metrics$rmse, 3), round(this_mod$metrics$r2, 3)))
  }
  else {
    NA_metrics <- tibble(ppt_id = ppt,
                         mae = NA,
                         rmse = NA,
                         r2 = NA)
    affect_metrics_intent <- affect_metrics_intent %>% add_row(NA_metrics)
  }
  
}

summary(affect_metrics_intent)
write.csv(affect_metrics_intent, "results/baseline_intent_affect.csv", row.names = FALSE)
summary(ar1_metrics_intent) # compare to AR1


