library(tidyr)
library(tidyverse)
library(dplyr) #this allows me to use the filter_all function
library(zoo) #this allows me to use the as.yearmon function
library(ggpubr)
library(ggplot2)
library(jtools)
library(lm.beta)
library(rstatix)

setwd("/Users/megsheppard/Desktop/OHBM R Analysis")

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




#rename the variables and combine to make overall grey matter volumes
#also whilst I'm doing this, times the volumes by the volume scaling factor
greyvol_insular <- (biobank$Volume_of_grey_matter_in_Insular_Cortex_.left..y+ biobank$Volume_of_grey_matter_in_Insular_Cortex_.right..y)*biobank$Volume_scaling_factor
greyvol_cingulate <- (biobank$Volume_of_grey_matter_in_Cingulate_Gyrus._anterior_division_.left..y+ biobank$Volume_of_grey_matter_in_Cingulate_Gyrus._anterior_division_.right..y)*biobank$Volume_scaling_factor
greyvol_parahippgyrus <- (biobank$Volume_of_grey_matter_in_Parahippocampal_Gyrus._anterior_division_.left..y + biobank$Volume_of_grey_matter_in_Parahippocampal_Gyrus._anterior_division_.right..y 
                          + biobank$Volume_of_grey_matter_in_Parahippocampal_Gyrus._posterior_division_.right..y + biobank$Volume_of_grey_matter_in_Parahippocampal_Gyrus._posterior_division_.left..y)*biobank$Volume_scaling_factor
greyvol_thal <- (biobank$Volume_of_grey_matter_in_Thalamus_.left..y+ biobank$Volume_of_grey_matter_in_Thalamus_.right..y)*biobank$Volume_scaling_factor
greyvol_hippocampus <- (biobank$Volume_of_grey_matter_in_Hippocampus_.left..y+biobank$Volume_of_grey_matter_in_Hippocampus_.right..y)*biobank$Volume_scaling_factor
greyvol_amygdala <- (biobank$Volume_of_grey_matter_in_Amygdala_.left..y+ biobank$Volume_of_grey_matter_in_Amygdala_.right..y)*biobank$Volume_scaling_factor
greyvol_striatum <- (biobank$Volume_of_grey_matter_in_Ventral_Striatum_.left..y + biobank$Volume_of_grey_matter_in_Ventral_Striatum_.right..y)*biobank$Volume_scaling_factor

#adds the new volume IDPS as columns in the biobank dataframe
biobank <- cbind(biobank, greyvol_insular)
biobank <- cbind(biobank, greyvol_cingulate)
biobank <- cbind(biobank, greyvol_parahippgyrus)
biobank <- cbind(biobank, greyvol_thal)
biobank <- cbind(biobank, greyvol_hippocampus)
biobank <- cbind(biobank, greyvol_striatum)
biobank <- cbind(biobank, greyvol_amygdala)

#sensitivity analyses
outliers.insula <- biobank %>% identify_outliers(greyvol_insular) %>% 
  filter(is.extreme==TRUE) %>% select(-c("is.outlier","is.extreme"))

biobank <- biobank %>% anti_join(outliers.insula)

outliers.cingulate <- biobank %>% identify_outliers(greyvol_cingulate) %>% 
  filter(is.extreme==TRUE) %>% select(-c("is.outlier","is.extreme"))

biobank <- biobank %>% anti_join(outliers.cingulate)

outliers.parahipp <- biobank %>% identify_outliers(greyvol_parahippgyrus) %>% 
  filter(is.extreme==TRUE) %>% select(-c("is.outlier","is.extreme"))

biobank <- biobank %>% anti_join(outliers.parahipp)

outliers.thal <- biobank %>% identify_outliers(greyvol_thal) %>% 
  filter(is.extreme==TRUE) %>% select(-c("is.outlier","is.extreme"))

biobank <- biobank %>% anti_join(outliers.thal)

outliers.hippocampus <- biobank %>% identify_outliers(greyvol_hippocampus) %>% 
  filter(is.extreme==TRUE) %>% select(-c("is.outlier","is.extreme"))

biobank <- biobank %>% anti_join(outliers.hippocampus)

outliers.amy <- biobank %>% identify_outliers(greyvol_amygdala) %>% 
  filter(is.extreme==TRUE) %>% select(-c("is.outlier","is.extreme"))

biobank <- biobank %>% anti_join(outliers.amy)

outliers.striatum <- biobank %>% identify_outliers(greyvol_striatum) %>% 
  filter(is.extreme==TRUE) %>% select(-c("is.outlier","is.extreme"))

biobank <- biobank %>% anti_join(outliers.striatum)

mean(biobank$child_stress)
sd(biobank$child_stress)
range(biobank$child_stress) 

#now doing multiple linear regression models that control for both age and sex 
#also running the associated confidence intervals
lm_amyg <- lm(greyvol_amygdala~biobank$child_stress + biobank$Age_at_visit + biobank$Sex_.0.F._1.M., data = biobank)
summary(lm_amyg)
confint(lm_amyg, level = 0.95)

summ(lm_amyg)
lm.beta(lm_amyg, complete.standardization = FALSE)

amy <- ggplot(data=biobank, aes(child_stress, greyvol_amygdala)) + geom_point()
thal <- ggplot(data=biobank, aes(child_stress, greyvol_thal)) + geom_point()
striatum <- ggplot(data=biobank, aes(child_stress, greyvol_striatum)) + geom_point()


lm_insular <- lm(greyvol_insular~biobank$child_stress + biobank$Age_at_visit + biobank$Sex_.0.F._1.M., data = biobank)
summary(lm_insular)
confint(lm_insular, level = 0.95)

lm_cingulate <- lm(greyvol_cingulate~biobank$child_stress + biobank$Age_at_visit + biobank$Sex_.0.F._1.M., data = biobank)
summary(lm_cingulate)
confint(lm_cingulate, level = 0.95)

lm_parahipp <- lm(greyvol_parahippgyrus~biobank$child_stress + biobank$Age_at_visit + biobank$Sex_.0.F._1.M., data = biobank)
summary(lm_parahipp)
confint(lm_parahipp, level = 0.95)

lm_thal <- lm(greyvol_thal~biobank$child_stress + biobank$Age_at_visit + biobank$Sex_.0.F._1.M., data = biobank)
summary(lm_thal)
confint(lm_thal, level = 0.95)

summ(lm_thal)
lm.beta(lm_thal, complete.standardization = FALSE)


lm_hippo <- lm(greyvol_hippocampus~biobank$child_stress + biobank$Age_at_visit + biobank$Sex_.0.F._1.M., data = biobank)
summary(lm_hippo) 
confint(lm_hippo, level = 0.95)

lm_striatum <- lm(greyvol_striatum~biobank$child_stress + biobank$Age_at_visit + biobank$Sex_.0.F._1.M., data = biobank)
summary(lm_striatum)
confint(lm_striatum, level = 0.95)

summ(lm_striatum)
lm.beta(lm_striatum, complete.standardization = FALSE)

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

#average and sd of CTQ score
av_CTQ <- mean(biobank$child_stress)
sd_CTQ <- sd(biobank$child_stress)

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

ggplot(biobank, aes(x=greyvol_amygdala, y=child_stress)) + geom_line()

library(rstatix)
outliers.amygdala<- biobank %>% identify_outliers(greyvol_amygdala) %>% 
  filter(is.extreme==TRUE) %>% select(-c("is.outlier","is.extreme"))

greyvol_amygdala <- biobank %>% anti_join(outliers.amygdala)

ggplot(greyvol_amygdala, aes(x=greyvol_amygdala, y=child_stress)) + geom_point() + geom_smooth()

outliers.striatum<- biobank %>% identify_outliers(greyvol_striatum) %>% 
filter(is.extreme==TRUE) %>% select(-c("is.outlier","is.extreme"))

greyvol_striatum <- biobank %>% anti_join(outliers.striatum)
ggplot(greyvol_striatum, aes(x=greyvol_striatum, y=child_stress)) + geom_point() + geom_smooth()
