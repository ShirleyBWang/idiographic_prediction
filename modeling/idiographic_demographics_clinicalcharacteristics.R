## Demographics and Clinical Characteristics for Idiographic Paper

##---- packages ----
library(tidyverse)
library(psych)
library(tidymodels)
library(tidyr)

##---- data -----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
u01_intent <- read.csv("../data/intent_affect_lagged_021623.csv")
u01_intent %>% group_by(ppt_id) %>% count() # 89 people 

sitbi_fc <- read.csv("../data/U01_Meditech_AllRC_Deidentified_1.24.2023.csv")
sitbi_mgh <- read.csv("../data/baseline_adult_cleaned_010923.csv")
demogs <- read.csv("../data/patient_demographics_010923.csv")

sitbi_fc_included <- sitbi_fc %>% filter(u01_ID %in% u01_intent$ppt_id)
sitbi_mgh_included <- sitbi_mgh %>% filter(ppt_id %in% u01_intent$ppt_id)
demogs_included <- demogs %>% filter(ppt_id %in% u01_intent$ppt_id)

##---- add column indicating yes/no included ----
sitbi_fc <- sitbi_fc %>% 
  mutate(idiographic_included = factor(ifelse(u01_ID %in% u01_intent$ppt_id, "Yes", "No")))

sitbi_mgh <- sitbi_mgh %>% 
  mutate(idiographic_included = factor(ifelse(ppt_id %in% u01_intent$ppt_id, "Yes", "No")))

demogs <- demogs %>% 
  mutate(idiographic_included = factor(ifelse(ppt_id %in% u01_intent$ppt_id, "Yes", "No")))

##---- SITBI recoding to merge ----
dim(sitbi_fc)
dim(sitbi_mgh)
names(sitbi_fc)
names(sitbi_mgh)

# suicidal thoughts: SITB_1 (si_yesno), SITB_5 (si_week), SITB_6 (si_month), SITB_7 (si_yearwks)
# suicide plan: SITB_10 (splan_yesno), SITB_14 (splan_week), SITB_15 (splan_month), SITB_16 (splan_yearwks)
# suicide attempt: SITB_19 (sa_yesno), SITB_23(sa_week), SITB_24 (sa_month), SITB_25 (sa_year)

sitbi_fc <- sitbi_fc %>% 
  select(u01_ID, idiographic_included, si_yesno, si_week, si_month, si_yearwks,
         splan_yesno, splan_week, splan_month, splan_yearwks,
         sa_yesno, sa_week, sa_month, sa_year) %>%
  rename(ppt_id = u01_ID)
sitbi_fc$sample <- "FC"

# do a bunch of renaming/recoding on the adult dataset
sitbi_mgh <- sitbi_mgh %>% 
  select(ppt_id, idiographic_included, SITB_1, SITB_5, SITB_6, SITB_7, 
         SITB_10, SITB_14, SITB_15, SITB_16,
         SITB_19, SITB_23, SITB_24, SITB_25) %>% 
  rename(si_yesno = SITB_1, si_week = SITB_5, si_month = SITB_6, si_yearwks = SITB_7, 
         splan_yesno = SITB_10, splan_week = SITB_14, splan_month = SITB_15, splan_yearwks = SITB_16,
         sa_yesno = SITB_19, sa_week = SITB_23, sa_month = SITB_24, sa_year = SITB_25)
sitbi_mgh$sample <- "MGH"
sitbi_mgh_recoded <- sitbi_mgh
sitbi_mgh_recoded$si_yesno <- recode(sitbi_mgh_recoded$si_yesno, "No" = 0, "Yes" = 1)
sitbi_mgh_recoded$si_week <- recode(sitbi_mgh_recoded$si_week, "0 days" = 0, "1 day" = 1, "2 days" = 2, 
                                    "3 days" = 3, "4 days" = 4, "5 days" = 5, "6 days" = 6, "7 days" = 7)

## recode specific instances
sitbi_mgh_recoded$si_month <- ifelse(sitbi_mgh_recoded$si_month == "25-30", 27, as.numeric(sitbi_mgh_recoded$si_month))
sitbi_mgh_recoded$si_yearwks <- ifelse(sitbi_mgh_recoded$si_yearwks == "10 or more", 10, as.numeric(sitbi_mgh_recoded$si_yearwks))

## no more specific instances
sitbi_mgh_recoded$splan_yesno <- recode(sitbi_mgh_recoded$splan_yesno, "No" = 0, "Yes" = 1)
sitbi_mgh_recoded$splan_week <- recode(sitbi_mgh_recoded$splan_week,"0 days" = 0, "1 day" = 1, 
                                       "2 days" = 2, "3 days" = 3, "4 days" = 4, "5 days" = 5, "6 days" = 6, "7 days" = 7)
sitbi_mgh_recoded$splan_month <- as.numeric(sitbi_mgh_recoded$splan_month)
sitbi_mgh_recoded$splan_yearwks <- as.numeric(sitbi_mgh_recoded$splan_yearwks)
sitbi_mgh_recoded$sa_yesno <- recode(sitbi_mgh_recoded$sa_yesno, "No" = 0, "Yes" = 1)

## merge the SITBI datasets 
sitbi <- sitbi_fc %>% add_row(sitbi_mgh_recoded)

##---- SITBI analyses ----
sitbi %>% 
  filter(idiographic_included == "Yes") %>% 
  summary()

# compare between groups
par(mfrow = c(3, 3))
plot(si_yearwks ~ idiographic_included, data = sitbi)
plot(si_month ~ idiographic_included, data = sitbi)
plot(si_week ~ idiographic_included, data = sitbi)

plot(splan_yearwks ~ idiographic_included, data = sitbi)
plot(splan_month ~ idiographic_included, data = sitbi)
plot(splan_week ~ idiographic_included, data = sitbi)

plot(sa_year ~ idiographic_included, data = sitbi)
plot(sa_month ~ idiographic_included, data = sitbi)
plot(sa_week ~ idiographic_included, data = sitbi)
par(mfrow = c(1, 1))

describeBy(sitbi, group = sitbi$idiographic_included)

t_test_results <- sitbi %>% 
  select(si_week, si_month, si_yearwks,
         splan_week, splan_month, splan_yearwks,
         sa_week, sa_month, sa_year) %>% 
  map(~ t.test(.x ~ sitbi$idiographic_included))
t_test_results

wilcox_test_results <- sitbi %>% 
  select(si_week, si_month, si_yearwks,
         splan_week, splan_month, splan_yearwks,
         sa_week, sa_month, sa_year) %>% 
  map(~ wilcox.test(.x ~ sitbi$idiographic_included))
wilcox_test_results

##---- demographic analyses ----
demogs_included %>% describe()
summary(demogs_included$Race_White)
summary(demogs_included$Race_Black)
summary(demogs_included$Race_Asian)
summary(as.factor(demogs_included$Hispanic))
summary(as.factor(demogs_included$Race))

## recode race to add AIAN/NHPI
demogs_included$Race_AIAN <- str_detect(demogs_included$Race, "Alaskan")
demogs_included$Race_NHPI <- str_detect(demogs_included$Race, "Hawaiian")

## recode gender in demogs
demogs <- demogs %>%
  mutate(GenderID = factor(GenderID))

demogs$GenderID <- demogs$GenderID %>% 
  fct_collapse("Non-binary/Gender non-conforming" = c("Non-binary/Gender non-conforming",
                                                      "Non-binary/ Gender-nonconforming"),
               "Other" = c("Other", "Other (please use text box to specify)"))

demogs$GenderID %>% summary()

## recode gender in demogs_included
demogs_included <- demogs_included %>%
  mutate(GenderID = factor(GenderID))

demogs_included$GenderID <- demogs_included$GenderID %>% 
  fct_collapse("Non-binary/Gender non-conforming" = c("Non-binary/Gender non-conforming",
                                                      "Non-binary/ Gender-nonconforming"),
               "Other" = c("Other", "Other (please use text box to specify)"))

demogs_included$GenderID %>% summary()

# compare results between groups
table(demogs$Race_White, demogs$idiographic_included) %>% chisq.test()
table(demogs$Race_Black, demogs$idiographic_included) %>% chisq.test()
table(demogs$Race_Asian, demogs$idiographic_included) %>% chisq.test()
table(demogs$Hispanic_Yes, demogs$idiographic_included) %>% chisq.test()
table(demogs$GenderID, demogs$idiographic_included) %>% chisq.test()

# write to csv
write.csv(demogs_included, "../data/demogs_included.csv", row.names = FALSE)
sitbi_included <- sitbi %>% 
  filter(idiographic_included == "Yes")
write.csv(sitbi_included, "../data/sitbi_included.csv", row.names = FALSE)


