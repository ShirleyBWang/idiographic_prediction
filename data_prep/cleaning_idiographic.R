## Data Cleaning for GP Models ##

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
conflict_prefer("lag", "dplyr")
`%notin%` <- Negate(`%in%`)

##---- load data ----
setwd("/Users/shirleywang/Documents/Papers/Matt/EMA/U01_Idiographic/u01_idiographic_analyses")
fulldat <- read.csv('data/ema_full_cleaned_110122.csv')
fulldat$phase %>% factor %>% summary()

##---- filter out daily surveys ----
dat <- fulldat %>% filter(survey_name != "Daily Survey") %>% 
  filter(day_in_study <= 84)
dat$phase %>% factor %>% summary()

dat %>%
  group_by(ppt_id) %>%
  summarise(days_in_study = max(day_in_study)) %>% View()

##---- time ----

## parse date and time
dat$start_time <- ymd_hms(dat$started_at)

## order by date/time within each sub
dat <- dat %>% arrange(ppt_id, start_time)

## add response number WITHIN each day and TOTAL number of responses
dat$resp_num_day <- rep(0, nrow(dat))
counter <- 1
for (i in unique(dat$ppt_id)){
  subdat <- subset(dat, ppt_id == i)
  for (j in unique(subdat$day_in_study)) {
    npings <- sum(subdat$day_in_study == j)
    dat$resp_num_day[counter:(counter + npings - 1)] <- 1:npings
    counter <- counter + npings
  }
}
summary(dat$resp_num_day)
dat %>% filter(resp_num_day > 6) %>% select(ppt_id) %>% unique() %>% nrow() # 148 ppt with > 6 EMA per day

dat <- dat %>% group_by(ppt_id) %>% mutate(resp_num_tot = n())

# dat %>% select(ppt_id, start_time, day_in_study, resp_num_day) %>% View() # seems to work! 

##---- clean: add time since entered study ----
dat_time <- dat

# parse date and time of ENTRANCE into study
dat_time$data_collection_start_parsed <- ymd_hms(dat_time$data_collection_start)

# get difftime between study entrance and present survey
dat_time$time_in_study <- difftime(dat_time$start_time, dat_time$data_collection_start_parsed, units = "hours")

## counts of missing urge and intent ratings
dat_time %>% select(sitb_si_intent, sitb_si_urge) %>% is.na() %>% summary() # approx 4500 missing, 66000 not

##---- create dataset for urge analyses ----
# at least 50 responses 
dat_urge <- dat_time %>% 
  group_by(ppt_id) %>%
  filter(sum(!is.na(sitb_si_urge)) >= 100) 

# let's start with urge only (first remove NAs, then re-create lagged and lead variables)
u01_urge <- dat_urge %>% filter(!is.na(sitb_si_urge))
u01_urge <- u01_urge %>% 
  group_by(ppt_id) %>% 
  mutate(time_lead = as.numeric(difftime(lead(start_time), start_time, units = "hours")),
         intent_lead = lead(sitb_si_intent),
         urge_lead = lead(sitb_si_urge), 
         lag_sisum = lag(sitb_si_sum), 
         lag_intent = lag(sitb_si_intent),
         lag_urge = lag(sitb_si_urge))

u01_urge %>% 
  group_by(ppt_id) %>% 
  summarise(time_count = sum(time_lead < 48, na.rm = TRUE)) %>% View()

urge_timelead_criteria <- u01_urge %>% 
  group_by(ppt_id) %>%
  summarise(time_lead_criteria = max(time_lead, na.rm = TRUE)) %>%
  filter(time_lead_criteria <= 72)

u01_urge <- u01_urge %>% 
  filter(ppt_id %in% urge_timelead_criteria$ppt_id)

lag_affect <- u01_urge %>% 
  group_by(ppt_id) %>% 
  select(all_of(affect_vars)) %>% 
  mutate_all(lag) %>% 
  rename_with( ~ str_replace(., "affect_", "lag_"))

urge_affect_lagged <- cbind(u01_urge, lag_affect[-1])
urge_affect_lagged %>% count() # 84 people 

urge_affect_lagged %>% 
  group_by(ppt_id) %>% 
  summarise(max = max(time_lead, na.rm = TRUE)) %>% View()

write.csv(urge_affect_lagged, "data/urge_affect_lagged_021623.csv", row.names = FALSE)
which(unique(intent_affect_lagged$ppt_id) %notin% unique(urge_affect_lagged$ppt_id)) # u01-mgh259 in intent but not urge
unique(intent_affect_lagged$ppt_id)[85]

##---- create dataset for intent analyses ----
dat_intent <- dat_time %>% 
  group_by(ppt_id) %>%
  filter(sum(!is.na(sitb_si_intent)) >= 100)

# let's start with urge only (first remove NAs, then re-create lagged and lead variables)
u01_intent <- dat_intent %>% filter(!is.na(sitb_si_intent))
u01_intent <- u01_intent %>% 
  group_by(ppt_id) %>% 
  mutate(time_lead = as.numeric(difftime(lead(start_time), start_time, units = "hours")),
         intent_lead = lead(sitb_si_intent),
         urge_lead = lead(sitb_si_urge), 
         lag_sisum = lag(sitb_si_sum), 
         lag_intent = lag(sitb_si_intent),
         lag_urge = lag(sitb_si_urge))

intent_timelead_criteria <- u01_intent %>% 
  group_by(ppt_id) %>%
  summarise(time_lead_criteria = max(time_lead, na.rm = TRUE)) %>%
  filter(time_lead_criteria <= 72)

u01_intent <- u01_intent %>% 
  filter(ppt_id %in% intent_timelead_criteria$ppt_id)

lag_affect <- u01_intent %>% 
  group_by(ppt_id) %>% 
  select(all_of(affect_vars)) %>% 
  mutate_all(lag) %>% 
  rename_with( ~ str_replace(., "affect_", "lag_"))

intent_affect_lagged <- cbind(u01_intent, lag_affect[-1])
intent_affect_lagged %>% count() # 89 people 

intent_affect_lagged %>% 
  group_by(ppt_id) %>% 
  summarise(max = max(time_lead, na.rm = TRUE)) %>% View()

write.csv(intent_affect_lagged, "data/intent_affect_lagged_021623.csv", row.names = FALSE)



