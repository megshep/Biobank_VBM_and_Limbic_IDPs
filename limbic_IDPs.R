library(tidyr)
library(tidyverse)
library(dplyr) #this allows me to use the filter_all function
library(zoo) #this allows me to use the as.yearmon function
library(ggpubr)
library(ggplot2)
library(jtools)
library(lm.beta)
library(rstatix)
library(ggeffects)
library(dplyr)
library(gridExtra)

meg_data <- read.csv("/Users/megsheppard/Desktop/Desktop - Meg’s MacBook Pro/OHBM R Analysis/data.csv")

meg_data

#this creates a dataset of just the non-cancerous illness codes, I've done this so I can remove any of which may be my exclusion criteria before re-merging it all back to the full dataset
data_exc <- meg_data[,c(2,291:319 )]

#this excludes all of my predetermined non-cancerous illness exclusion criteria using the Biobank codes for each criterion.
data_exc_all <- data_exc %>% filter_all(any_vars(. %in% 
              c(1081, 1082, 1083, 1086, 1240, 1245, 1246, 
                1247, 1249, 1250, 1258, 1259,1261, 1262, 
                1263,1264, 1266,1397, 1434, 1491, 1495,
                1524, 1583, 1659,9999)))


#this creates a new dataframe which has excluded all the individuals with any of the non-cancerous illness codes above (excludes NCI dataframe from medical dataframe)
biobank_exc <- meg_data[!(meg_data$PID %in% data_exc_all$PID),] #this matched people on eID as it will only merge in participants present in both dataframes so will isolate and remove excluded individuals

#Getting rid of participants who don't have answers to the CTQ
#creating a dataset containing only pid and CTQ scores]
CTQ_exc <- biobank_exc[,c(2, 13:17)]

#deleting rows where they are all NA - aka participants don't have any answers
CTQ_exc <- CTQ_exc[!apply(CTQ_exc[2:ncol(CTQ_exc)], MARGIN =1, 
                          FUN = function(x) {all(is.na(x))}),]

#add back into main data set - merging by PID' ensures that only rows where pids appear in both datasets are added back in(but this has created a new column without removing the first). 
biobank_CTQexc <- merge(biobank_exc, CTQ_exc, by = "PID")

#need to now repeat this for the IDPs to get rid of anyone who doesnt have a full set of IDPS 
#creating a variable containing only eid and idps 
vol_ex <- biobank_CTQexc[,c(1, 272:287)]

#deleting rows where they are all NA - aka participants don't have any IDPs
vol_ex <- vol_ex[!apply(vol_ex[2:ncol(vol_ex)], MARGIN =1, 
                        FUN = function(x) {all(is.na(x))}),]

#add back into main data set - merging by 'eid' ensures that only rows where eids appear in both datasets are added back in. 
#ukb_excluded <- merge(main_dataset, img_ex, by = "eid")
biobank_excluded <- merge(biobank_CTQexc, vol_ex, by = "PID")

#this deletes/removes the original columns from the dataframe 
biobank_excluded <- biobank_excluded[, -c(13:17, 272:287)]

#renaming the CTQ variables from the dataframe 
trauma_1 <- biobank_excluded$Traumatic_child_events_1.y
trauma_2 <- biobank_excluded$Traumatic_child_events_2.y
trauma_3 <- biobank_excluded$Traumatic_child_events_3.y
trauma_4 <- biobank_excluded$Traumatic_child_events_4.y
trauma_5 <- biobank_excluded$Traumatic_child_events_5.y
  
#now we need to add up these scores for each participant for their overall CTQ scores - need to add to dataframe 
CTQ <- trauma_1 + trauma_2 + trauma_3 + trauma_4 + trauma_5

#this adds the CTQ as a column in my new and hopefully final dataframe
biobank <- cbind(biobank_excluded, CTQ)


#rename the variable and combine to make overall grey matter volume
#also whilst I'm doing this, times the volume by the volume scaling factor
greyvol_amygdala <- (biobank$Volume_of_grey_matter_in_Amygdala_.left..y+ biobank$Volume_of_grey_matter_in_Amygdala_.right..y)*biobank$Volume_scaling_factor

#adds the new volume IDP as a column in the biobank dataframe
biobank <- cbind(biobank, greyvol_amygdala)

#sensitivity analysis

outliers.amy <- biobank %>% identify_outliers(greyvol_amygdala) %>% 
  filter(is.extreme==TRUE) %>% select(-c("is.outlier","is.extreme"))

biobank <- biobank %>% anti_join(outliers.amy)

#work out the average, sd and range of CTQ scores within my sample
mean(biobank$child_stress)
sd(biobank$child_stress)
range(biobank$child_stress) 


##############demographics

#working out the average age and sd of the whole sample
av_age_all <- mean(biobank$Age_at_visit)
sd_age_all <- sd(biobank$Age_at_visit)

#working out the average age and sd of just the women within the sample
female <- biobank %>% filter (Sex_.0.F._1.M. == 0) %>%
  summarise(av_age = mean(Age_at_visit))

female_sd <- biobank %>% filter (Sex_.0.F._1.M. == 0) %>%
  summarise(sd_age = sd(Age_at_visit))

#working out the average age and sd of just the men within the sample
male <- biobank %>% filter (Sex_.0.F._1.M. == 1) %>%
  summarise(av_age = mean (Age_at_visit))

male_sd <- biobank %>% filter (Sex_.0.F._1.M. == 1) %>%
  summarise(sd_age = sd(Age_at_visit))

#number of ethnicity entries
ethnicity_1001 <- biobank %>% count(ethnicity == 1001)
ethnicity_1002 <- biobank %>% count(ethnicity == 1002)
ethnicity_1003 <- biobank %>% count(ethnicity == 1003)
ethnicity_4001 <- biobank %>% count(ethnicity == 4001)
ethnicity_4002 <- biobank %>% count(ethnicity == 4002)
ethnicity_4003 <- biobank %>% count(ethnicity == 4003)
ethnicity_5 <- biobank %>% count(ethnicity == 5)
ethnicity_1 <- biobank %>% count(ethnicity == 1)
ethnicity_2 <- biobank %>% count(ethnicity == 2)
ethnicity_4 <- biobank %>% count(ethnicity == 4)
ethnicity_6 <- biobank %>% count(ethnicity == 6)
ethnicity_2003 <- biobank %>% count(ethnicity == 2003)
ethnicity_2001 <- biobank %>% count(ethnicity == 2001)
ethnicity_2002 <- biobank %>% count(ethnicity == 2002)
ethnicity_3004 <- biobank %>% count(ethnicity == 3004)
ethnicity_3001 <- biobank %>% count(ethnicity == 3001)
ethnicity_3002 <- biobank %>% count(ethnicity == 3002)
ethnicity_3003<- biobank %>% count(ethnicity == 3003)
ethnicity_2004 <- biobank %>% count(ethnicity == 2004)
ethnicity_neg3 <- biobank %>% count(ethnicity == -3)
ethnicity_neg1 <- biobank %>% count(ethnicity == -1)

biobank$Age_c <- biobank$Age_at_visit - mean(biobank$Age_at_visit, na.rm = TRUE)

# Modelling

#initial regression model to look at the associations between ELS and amygdala volume ±95%CIs
lm_amyg <- lm(greyvol_amygdala~biobank$child_stress + biobank$Age_at_visit + biobank$Sex_.0.F._1.M., data = biobank)
summary(lm_amyg)
confint(lm_amyg, level = 0.95)

summ(lm_amyg)
lm.beta(lm_amyg, complete.standardization = FALSE)

summary(lm_amyg)

#this is the centered age model looking at the interaction of cta and age polynomially and linearly
#this is now an exploratory trajectory analysis (out of interest, limited due to cross-sectional nature)

#defining the interaction term to look at the effect of ELS on trajectory within this sample
lm_amyg <- lm( greyvol_amygdala ~ child_stress * Age_c + Sex_.0.F._1.M., data = biobank)

#defining the interaction term polynomially 
lm_amyg_nonlinear <- lm( greyvol_amygdala ~ child_stress * (Age_c + I(Age_c^2)) + Sex_.0.F._1.M., data = biobank)
summary(lm_amyg_nonlinear)

#visualising the differences (or lack there of) between polynomial trajectories of 2±sd from the mean
pred <- ggpredict(
  lm_amyg_nonlinear,
  terms = c("Age_c [all]", "child_stress [-1,4]")
)

ggplot(pred,
       aes(x = x, y = predicted,
           colour = group,
           fill = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              alpha = .2,
              colour = NA) +
  labs(
    x = "Centered age",
    y = "Predicted amygdala volume",
    colour = "CTQ",
    fill = "CTQ"
  ) +
  theme_classic()


### now looking at age bins to see if there are any age specific effects that are lost in the wider model
biobank$age_group <- cut(biobank$Age_at_visit,breaks = c(49, 55, 60, 65, 70), include.lowest = TRUE)

#running a model for each of the bins,
lm_bins <- lm( greyvol_amygdala ~ child_stress * age_group + Sex_.0.F._1.M., data = bioban)
summary(lm_bins)

#generate predicted effects inc confidence intervals that can be read into the plot
eff <- ggpredict(lm_amyg, terms = "child_stress")

#this plots the overall regression line and includes a ribbon for the 95% confidence intervals from the predicted data defined above
p1 <- ggplot() + geom_line(data = eff, aes(x = x, y = predicted), linewidth = 1) +
  geom_ribbon(data = eff, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  labs(x = "Childhood stress", y = "Amygdala grey matter volume") +
  theme_classic(base_size = 12) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

#plotting a standardised beta plot to outline the contributions of each of the variables within th model
coefs <- summary(lm_amyg)$coefficients
coefs <- as.data.frame(coefs)
coefs$term <- rownames(coefs)

# remove intercept
coefs <- coefs[coefs$term != "(Intercept)", ]

# standardise using sds to help create partially standardized coefficient to allow us to cmpare the predictors
sd_y <- sd(biobank$greyvol_amygdala, na.rm = TRUE)
sd_x_child <- sd(biobank$child_stress, na.rm = TRUE)
sd_x_age   <- sd(biobank$Age_at_visit, na.rm = TRUE)
sd_x_sex   <- sd(biobank$Sex_.0.F._1.M., na.rm = TRUE)

#puts the standardised betas and standard deviations into a table to plot 
coefs <- coefs %>%  mutate( beta_std = Estimate * c(sd_x_child, sd_x_age, sd_x_sex), se_std   = `Std. Error` * c(sd_x_child, sd_x_age, sd_x_sex),
    term = case_when( term == "child_stress" ~ "Childhood Stress", term == "Age_at_visit" ~ "Age", term == "Sex_.0.F._1.M." ~ "Sex",
      TRUE ~ term))

#plot the information ±95%CIs
p2 <- ggplot(coefs, aes(x = reorder(term, beta_std), y = beta_std)) +
  geom_point(size = 2) + geom_errorbar(aes( ymin = beta_std - 1.96 * se_std,ymax = beta_std + 1.96 * se_std),
    width = 0.2) +
  coord_flip() +
  labs(x = NULL, y = "Standardised beta (β, 95% CI)") + theme_classic(base_size = 12)

#combine both figures to one final overall figure.
final_fig <- grid.arrange(p1, p2, ncol = 2)

