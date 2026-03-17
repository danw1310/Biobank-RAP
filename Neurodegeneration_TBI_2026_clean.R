##########Install Packages##################################

install.packages("dplyr")
install.packages("tidyverse")
install.packages("survival")
install.packages("survminer")
install.packages("epiR")
install.packages("stringr")
install.packages("naniar")
install.packages("tableone")
install.packages("AF")
install.packages("cmprsk")
install.packages("gridExtra")
install.packages("cowplot")
install.packages("emmeans")
install.packages("rstatix")
install.packages("broom")
install.packages("sjPlot")
install.packages("sjmisc")
install.packages("sjlabelled")
install.packages("mice")
install.packages("miceadds")
install.packages("car")
library(car)
library(miceadds)
library(mice)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(rstatix)
library(broom)
library(emmeans)
library(cowplot)
library(gridExtra)
library(cmprsk)
library(tableone)
library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(epiR)
library(stringr)
library(naniar)
library(AF)


# Define the function to export Cox model summary to CSV
export_cox_summary_to_csv <- function(cox_model, file_path) {
  # Check if the input is a coxph model
  if (!inherits(cox_model, "coxph")) {
    stop("The input model is not a Cox proportional hazards model (coxph).")
  }
  
  # Get the summary of the coxph model
  summary_cox_model <- summary(cox_model)
  
  # Extract relevant information from the summary
  coefficients <- summary_cox_model$coefficients
  conf.int <- summary_cox_model$conf.int
  covariate_names <- rownames(coefficients)
  
  # Combine covariate names, hazard ratios with confidence intervals, and p-values into a data frame
  result <- data.frame(
    covariate = covariate_names,
    hazard_ratio = sprintf("%.2f (%.2f, %.2f)", 
                           coefficients[, "exp(coef)"], 
                           conf.int[, "lower .95"], 
                           conf.int[, "upper .95"]),
    p_value = coefficients[, "Pr(>|z|)"]
  )
  
  # Write the result to a CSV file
  write.csv(result, file = file_path, row.names = FALSE)
  
  # Print a message indicating where the file has been saved
  print(paste("CSV file has been created at:", file_path))
}


#############Dataframes######################################

Head_injured <- read.csv("~/Daniel_Whitehouse/Head_injured_cohort/output/complete_head_injured_Core.csv",na.strings = c("", "NA"))
TBI <- read.csv("~/Daniel_Whitehouse/Head_injured_cohort/output/complete_head_injured_Core.csv",na.strings = c("", "NA"))
Nero <- read.csv("~/Daniel_Whitehouse/Neurodegeneration/input/Neurodegeneration.csv",na.strings = c("", "NA"))
Age <- read.csv("~/Daniel_Whitehouse/Head_injured_cohort/input/Age.csv",na.strings = c("", "NA"))
Demographics <- read.csv("~/Daniel_Whitehouse/Neurodegeneration/input/Demographics_neurodegeneration.csv",na.strings = c("", "NA"))
recruitment <- read.csv("~/Daniel_Whitehouse/Neurodegeneration/input/Recruitment_timings.csv",na.strings = c("", "NA"))
cholesterol <- read.csv("~/Daniel_Whitehouse/Neurodegeneration/input/Cholesterol_participant.csv",na.strings = c("", "NA"))

withdrawn <- read.csv("~/Daniel_Whitehouse/Head_injured_cohort/output/Withdrawn_update.csv", na.strings = "")

Nero <- subset(Nero, !Nero$eid%in%withdrawn$eid)

##Co-morbidity

Co_morbid <- read.csv("~/Daniel_Whitehouse/Neurodegeneration/input/Co-morbidities_DAG.csv",na.strings = c("", "NA"))
Death_reg <- read.csv("~/Daniel_Whitehouse/Neurodegeneration/input/Death_reg.csv",na.strings = c("", "NA"))
Genetics <- read.csv("~/Daniel_Whitehouse/Neurodegeneration/input/PRS_dementia_PD.csv",na.strings = c("", "NA"))

Non_withdrawn_participants <- Nero %>% select(eid)

Head_injured <- Head_injured %>% filter(eid %in% Non_withdrawn_participants$eid)
#TBI <- TBI %>% filter(eid %in% Non_withdrawn_participants$eid)
Nero <- Nero %>% filter(eid %in% Non_withdrawn_participants$eid)
Age <- Age %>% filter(eid %in% Non_withdrawn_participants$eid)
Demographics <- Demographics %>% filter(eid %in% Non_withdrawn_participants$eid)
recruitment <- recruitment %>% filter(eid %in% Non_withdrawn_participants$eid)
Co_morbid <- Co_morbid %>% filter(eid %in% Non_withdrawn_participants$eid)
Death_reg <- Death_reg %>% filter(eid %in% Non_withdrawn_participants$eid)
Genetics <- Genetics %>% filter(eid %in% Non_withdrawn_participants$eid)
cholesterol <- cholesterol %>% filter(eid %in% Non_withdrawn_participants$eid)
###Demographics clean

Demographics <- merge(Genetics, Demographics, by="eid")

key <- c(
  "African" = "Black British",
  "Any other Asian background" = "Asian British",
  "Any other Black background" = "Black British",
  "Any other mixed background" = "Other",
  "Any other white background" = "White British",
  "Asian or Asian British" = "Asian British",
  "Bangladeshi" = "Asian British",
  "Black or Black British" = "Black British",
  "British" = "White British",
  "Caribbean" = "Black British",
  "Chinese" = "Asian British",
  "Do not know" = "Other",
  "Indian" = "Asian British",
  "Irish" = "White British",
  "Mixed" = "Other",
  "Other ethnic group" = "Other",
  "Pakistani" = "Asian British",
  "Prefer not to answer" = "Other",
  "White" = "White British",
  "White and Asian" = "Asian British",
  "White and Black African" = "Black British",
  "White and Black Caribbean" = "Black British"
)

# Update values in the column based on the key
Demographics$p21000_i0 <- key[Demographics$p21000_i0]

table(Demographics$p21000_i0)

table(Demographics$p40022)
key_location <- c(
  "Hospital Episode Statistics from England" = "Hospital Episode Statistics from England",
  "Hospital Episode Statistics from England|Patient Episode Database for Wales"="Patient Episode Database for Wales",
  "Hospital Episode Statistics from England|Patient Episode Database for Wales|Scottish Morbidity Records"= "Scottish Morbidity Records",
  "Hospital Episode Statistics from England|Scottish Morbidity Records"="Scottish Morbidity Records",
  "Patient Episode Database for Wales"="Patient Episode Database for Wales",
  "Patient Episode Database for Wales|Scottish Morbidity Records" = "Scottish Morbidity Records",
  "Scottish Morbidity Records" = "Scottish Morbidity Records")
Demographics$p40022 <- key_location[Demographics$p40022]
table(Demographics$p40022)

Demographics$obesity <- ifelse(Demographics$p21001_i0 >= 30, 1,0)

Education <- Demographics %>% select(c(eid, p6138_i0))
Education_wide <- separate_wider_delim(data = Education, cols = p6138_i0, delim = "|", names_sep = "_", too_few = "align_start")

#College or University degree = college
#Other professional qualifications eg: nursing, teaching = college
#A levels/AS levels or equivalent = secondary
#O levels/GCSEs or equivalent = secondary
#CSEs or equivalent = secondary
#NVQ or HND or HNC or equivalent = college
#None of the above = primary
#Prefer not to answer

Education_wide$p6138_i0_1 <- str_replace(Education_wide$p6138_i0_1,'College or University degree','College')
Education_wide$p6138_i0_1 <- str_replace(Education_wide$p6138_i0_1,'Other professional qualifications eg: nursing, teaching','Training')
Education_wide$p6138_i0_1 <- str_replace(Education_wide$p6138_i0_1,'A levels/AS levels or equivalent','Secondary')
Education_wide$p6138_i0_1 <- str_replace(Education_wide$p6138_i0_1,'O levels/GCSEs or equivalent','Secondary')
Education_wide$p6138_i0_1 <- str_replace(Education_wide$p6138_i0_1,'CSEs or equivalent','Secondary')
Education_wide$p6138_i0_1 <- str_replace(Education_wide$p6138_i0_1,'NVQ or HND or HNC or equivalent','Training')
Education_wide$p6138_i0_1 <- str_replace(Education_wide$p6138_i0_1,'None of the above','Primary')

Education_wide$p6138_i0_2 <- str_replace(Education_wide$p6138_i0_2,'College or University degree','College')
Education_wide$p6138_i0_2 <- str_replace(Education_wide$p6138_i0_2,'Other professional qualifications eg: nursing, teaching','Training')
Education_wide$p6138_i0_2 <- str_replace(Education_wide$p6138_i0_2,'A levels/AS levels or equivalent','Secondary')
Education_wide$p6138_i0_2 <- str_replace(Education_wide$p6138_i0_2,'O levels/GCSEs or equivalent','Secondary')
Education_wide$p6138_i0_2 <- str_replace(Education_wide$p6138_i0_2,'CSEs or equivalent','Secondary')
Education_wide$p6138_i0_2 <- str_replace(Education_wide$p6138_i0_2,'NVQ or HND or HNC or equivalent','Training')
Education_wide$p6138_i0_2 <- str_replace(Education_wide$p6138_i0_2,'None of the above','Primary')

Education_wide$p6138_i0_3 <- str_replace(Education_wide$p6138_i0_3,'College or University degree','College')
Education_wide$p6138_i0_3 <- str_replace(Education_wide$p6138_i0_3,'Other professional qualifications eg: nursing, teaching','Training')
Education_wide$p6138_i0_3 <- str_replace(Education_wide$p6138_i0_3,'A levels/AS levels or equivalent','Secondary')
Education_wide$p6138_i0_3 <- str_replace(Education_wide$p6138_i0_3,'O levels/GCSEs or equivalent','Secondary')
Education_wide$p6138_i0_3 <- str_replace(Education_wide$p6138_i0_3,'CSEs or equivalent','Secondary')
Education_wide$p6138_i0_3 <- str_replace(Education_wide$p6138_i0_3,'NVQ or HND or HNC or equivalent','Training')
Education_wide$p6138_i0_3 <- str_replace(Education_wide$p6138_i0_3,'None of the above','Primary')

Education_wide$p6138_i0_4 <- str_replace(Education_wide$p6138_i0_4,'College or University degree','College')
Education_wide$p6138_i0_4 <- str_replace(Education_wide$p6138_i0_4,'Other professional qualifications eg: nursing, teaching','Training')
Education_wide$p6138_i0_4 <- str_replace(Education_wide$p6138_i0_4,'A levels/AS levels or equivalent','Secondary')
Education_wide$p6138_i0_4 <- str_replace(Education_wide$p6138_i0_4,'O levels/GCSEs or equivalent','Secondary')
Education_wide$p6138_i0_4 <- str_replace(Education_wide$p6138_i0_4,'CSEs or equivalent','Secondary')
Education_wide$p6138_i0_4 <- str_replace(Education_wide$p6138_i0_4,'NVQ or HND or HNC or equivalent','Training')
Education_wide$p6138_i0_4 <- str_replace(Education_wide$p6138_i0_4,'None of the above','Primary')

Education_wide$p6138_i0_5 <- str_replace(Education_wide$p6138_i0_5,'College or University degree','College')
Education_wide$p6138_i0_5 <- str_replace(Education_wide$p6138_i0_5,'Other professional qualifications eg: nursing, teaching','Training')
Education_wide$p6138_i0_5 <- str_replace(Education_wide$p6138_i0_5,'A levels/AS levels or equivalent','Secondary')
Education_wide$p6138_i0_5 <- str_replace(Education_wide$p6138_i0_5,'O levels/GCSEs or equivalent','Secondary')
Education_wide$p6138_i0_5 <- str_replace(Education_wide$p6138_i0_5,'CSEs or equivalent','Secondary')
Education_wide$p6138_i0_5 <- str_replace(Education_wide$p6138_i0_5,'NVQ or HND or HNC or equivalent','Training')
Education_wide$p6138_i0_5 <- str_replace(Education_wide$p6138_i0_5,'None of the above','Primary')

Education_wide$p6138_i0_6 <- str_replace(Education_wide$p6138_i0_6,'College or University degree','College')
Education_wide$p6138_i0_6 <- str_replace(Education_wide$p6138_i0_6,'Other professional qualifications eg: nursing, teaching','Training')
Education_wide$p6138_i0_6 <- str_replace(Education_wide$p6138_i0_6,'A levels/AS levels or equivalent','Secondary')
Education_wide$p6138_i0_6 <- str_replace(Education_wide$p6138_i0_6,'O levels/GCSEs or equivalent','Secondary')
Education_wide$p6138_i0_6 <- str_replace(Education_wide$p6138_i0_6,'CSEs or equivalent','Secondary')
Education_wide$p6138_i0_6 <- str_replace(Education_wide$p6138_i0_6,'NVQ or HND or HNC or equivalent','Training')
Education_wide$p6138_i0_6 <- str_replace(Education_wide$p6138_i0_6,'None of the above','Primary')

Education_wide$College <- ifelse((Education_wide$p6138_i0_1=="College" | 
                                    Education_wide$p6138_i0_2=="College"| 
                                    Education_wide$p6138_i0_3=="College"|
                                    Education_wide$p6138_i0_4=="College"|
                                    Education_wide$p6138_i0_5=="College"|
                                    Education_wide$p6138_i0_6=="College"),"College", "")

Education_wide$Secondary <- ifelse((Education_wide$p6138_i0_1=="Secondary" | 
                                      Education_wide$p6138_i0_2=="Secondary"| 
                                      Education_wide$p6138_i0_3=="Secondary"|
                                      Education_wide$p6138_i0_4=="Secondary"|
                                      Education_wide$p6138_i0_5=="Secondary"|
                                      Education_wide$p6138_i0_6=="Secondary"),"Secondary", "")

Education_wide$Primary <- ifelse((Education_wide$p6138_i0_1=="Primary" | 
                                    Education_wide$p6138_i0_2=="Primary"| 
                                    Education_wide$p6138_i0_3=="Primary"|
                                    Education_wide$p6138_i0_4=="Primary"|
                                    Education_wide$p6138_i0_5=="Primary"|
                                    Education_wide$p6138_i0_6=="Primary"),"Primary", "")

Education_wide$Training <- ifelse((Education_wide$p6138_i0_1=="Training" | 
                                     Education_wide$p6138_i0_2=="Training"| 
                                     Education_wide$p6138_i0_3=="Training"|
                                     Education_wide$p6138_i0_4=="Training"|
                                     Education_wide$p6138_i0_5=="Training"|
                                     Education_wide$p6138_i0_6=="Training"),"Training", "")


Education_wide$highest_level <- paste(Education_wide$College, Education_wide$Secondary, Education_wide$Primary, Education_wide$Training)
Education_wide$highest_level <- ifelse(grepl("College",Education_wide$highest_level),'College', Education_wide$highest_level)    
Education_wide$highest_level <- ifelse(grepl("Training",Education_wide$highest_level),'Training', Education_wide$highest_level)
Education_wide$highest_level <- ifelse(grepl("Secondary",Education_wide$highest_level),'Secondary', Education_wide$highest_level) 
Education_wide$highest_level <- ifelse(grepl("Primary",Education_wide$highest_level),'Primary', Education_wide$highest_level) 
Education_wide$highest_level <- ifelse(grepl("NA",Education_wide$highest_level),NA, Education_wide$highest_level)

Highest_education <- Education_wide %>% select(eid, highest_level)

Demographics <- merge(Demographics, Highest_education, by="eid", all.x = T)

Demographics$Lives_alone <- ifelse(Demographics$p709_i0 == "1", "Lives Alone",
                                   ifelse(Demographics$p709_i0 == "Do not know", "Do not know",
                                          ifelse(Demographics$p709_i0 == "Prefer not to answer", "Prefer not to answer",
                                                 ifelse(is.na(Demographics$p709_i0), NA, "Lives with others"))))                                        

Demographics$Lives_alone <- ifelse(Demographics$p670_i0 == "Sheltered accommodation" | Demographics$p670_i0 == "Care home", 
                                   "Lives with others", 
                                   Demographics$Lives_alone)

Demographics <- Demographics %>%
  mutate_all(~replace(., . %in% c("Do not know", "Prefer not to answer"), NA))

Demographics$smoking <- ifelse(Demographics$p20116_i0=="Previous" & Demographics$p1239_i0=="No", "Previous", Demographics$p1239_i0)
Demographics$smoking[Demographics$smoking == "NA"] <- NA
Demographics$smoking <- factor(Demographics$smoking,
                               levels = c("No", "Only occasionally", "Previous", "Yes, on most or all days"),
                               labels = c("Never", "Occasional", "Previous", "Frequent"))


Demographics$alcohol <- ifelse(Demographics$p20117_i0=="Previous" & Demographics$p1558_i0=="Never", "Previous", Demographics$p1558_i0)
Demographics$alcohol[Demographics$alcohol == "NA"] <- NA
Demographics$alcohol <- factor(Demographics$alcohol,
                               levels = c("Never", "Previous", "Special occasions only", "One to three times a month", "Once or twice a week", "Three or four times a week", "Daily or almost daily"),
                               labels = c("Never", "Previous", "Light", "Light", "Light", "Moderate", "Heavy" ))

#Age at recruitment
Demographics <- merge(Demographics, Age, by="eid", all.x = T)

Demographics$p53_i0<-as.factor(Demographics$p53_i0)
Demographics$p53_i0<-strptime(Demographics$p53_i0,format="%Y-%m-%d")
Demographics$p53_i0<-as.Date(Demographics$p53_i0,format="%Y-%m-%d")

Demographics$Date_of_Birth_imputed<-as.factor(Demographics$Date_of_Birth_imputed)
Demographics$Date_of_Birth_imputed<-strptime(Demographics$Date_of_Birth_imputed,format="%Y-%m-%d")
Demographics$Date_of_Birth_imputed<-as.Date(Demographics$Date_of_Birth_imputed,format="%Y-%m-%d")

Demographics$Age_at_recruitment <- difftime(Demographics$p53_i0,Demographics$Date_of_Birth_imputed,units=c("weeks"))
Demographics$Age_at_recruitment <-Demographics$Age_at_recruitment/52.1429
Demographics$Age_at_recruitment<-gsub("weeks","",as.character(Demographics$Age_at_recruitment))
Demographics$Age_at_recruitment<-as.numeric(Demographics$Age_at_recruitment)

Age_recruitment_df <- Demographics %>% select(eid, Age_at_recruitment)

cholesterol$p30690_i0 <- as.numeric(cholesterol$p30690_i0)
cholesterol$p30780_i0 <- as.numeric(cholesterol$p30780_i0)
cholesterol$raised_lipids <- ifelse(cholesterol$p30780_i0>=4.1,1,0)
cholesterol$combine_self <- ifelse(is.na(cholesterol$p6177_i0), cholesterol$p6153_i0, cholesterol$p6177_i0)

cholesterol <- cholesterol %>%
  mutate_all(~replace(., . %in% c("Do not know", "Prefer not to answer"), NA))

check_hypertension_medication <- function(text) {
  if (is.na(text)) {
    return(NA)
  } else if (grepl("Blood pressure medication", text, ignore.case = T)) {
    return(1)
  } else {
    return(0)
  }
}

check_cholesterol_medicaiton <- function(text) {
  if (is.na(text)) {
    return(NA)
  } else if (grepl("Cholesterol lowering medication", text, ignore.case = T)) {
    return(1)
  } else {
    return(0)
  }
}

check_insulin_medication <- function(text) {
  if (is.na(text)) {
    return(NA)
  } else if (grepl("Insulin", text, ignore.case = T)) {
    return(1)
  } else {
    return(0)
  }
}


cholesterol$Hypertension_meds <- sapply(cholesterol$combine_self, check_hypertension_medication)
cholesterol$Cholesterol_meds<- sapply(cholesterol$combine_self, check_cholesterol_medicaiton)
cholesterol$Insulin_meds<- sapply(cholesterol$combine_self, check_insulin_medication)

Co_morbid <- merge(Co_morbid, cholesterol, by="eid", all.x = T)


###Family history

Family_history <- read.csv("~/Daniel_Whitehouse/Neurodegeneration/input/Family_history_participant.csv",na.strings = c("", "NA"))

Family_history$DemFH <- ifelse(grepl("Alzheimer's disease/dementia", Family_history$p20107_i0, fixed = TRUE) | 
                                 grepl("Alzheimer's disease/dementia", Family_history$p20110_i0, fixed = TRUE), 1, 0)

Family_history$PDFH <- ifelse(grepl("Parkinson's disease", Family_history$p20107_i0, fixed = TRUE) | 
                                grepl("Parkinson's disease", Family_history$p20110_i0, fixed = TRUE), 1, 0)

Demographics <- merge(Demographics, Family_history, by="eid", all.x = T)

##Genetics

APOE <- read.csv("~/Daniel_Whitehouse/Neurodegeneration/input/APOE_status.csv")
APOE <- subset(APOE, !APOE$APOE_Genotype=="ε2ε4")
Demographics <- merge(Demographics, APOE, by.x="eid", by.y = "FID", all.x = T)

#Audio

Audio_Visual <- read.csv("~/Daniel_Whitehouse/Neurodegeneration/input/Vision_Audio_participant.csv",na.strings = c("", "NA"))
Audio_Visual <- Audio_Visual %>%
  mutate_all(~replace(., . %in% c("Prefer not to answer"), NA))

Audio_Visual$hearing_loss <- ifelse(Audio_Visual$p2247_i0=="Yes"| Audio_Visual$p2247_i0=="I am completely deaf", 1, 0)
Audio_Visual$visual_loss <- ifelse(Audio_Visual$p2207_i0=="Yes", 1, 0)

Demographics <- merge(Demographics, Audio_Visual, by="eid", all.x = T)

##Co-morbiditiy cleaning 

Co_morbid <- merge(Co_morbid, Age, by="eid", all.x = T)
Co_morbid <- merge(Co_morbid, Age_recruitment_df, by="eid", all.x = T)
Co_morbid$Age_at_recruitment <- as.numeric(Co_morbid$Age_at_recruitment)
Co_morbid$Date_of_Birth_imputed<-as.factor(Co_morbid$Date_of_Birth_imputed)
Co_morbid$Date_of_Birth_imputed<-strptime(Co_morbid$Date_of_Birth_imputed,format="%Y-%m-%d")
Co_morbid$Date_of_Birth_imputed<-as.Date(Co_morbid$Date_of_Birth_imputed,format="%Y-%m-%d")

#Hypertension: p131286,p131288,p131290,p131292,p131294
#Diabetes: p130706,p130708,p130710,p130712,p130714
#Cerebrovascular: p42006, p131360,p131362,p131364,p131366, p131368,p131374,p131376, p131378
#Cardiovascular: p42000, p131296,p131298,p131300,p131302,p131304,p131306
#Depression: p130894, p130896,p130898, p130900, p130902

date_columns <- c("p131286","p131288","p131290","p131292",
                  "p131294","p130706","p130708","p130710","p130712","p130714", 
                  "p42006", "p131360","p131362","p131364","p131366", "p131368","p131374","p131376", "p131378", 
                  "p42000", "p131296","p131298","p131300","p131302","p131304","p131306",
                  "p130894", "p130896","p130898", "p130900", "p130902"
)

for (col in date_columns) {
  
  Co_morbid[[col]] <- as.factor(Co_morbid[[col]])
  
  Co_morbid[[col]] <- strptime(Co_morbid[[col]], format="%Y-%m-%d")
  
  Co_morbid[[col]] <- as.Date(Co_morbid[[col]], format="%Y-%m-%d")
}

##Hyperlipidemia - difined by either a total cholesterol >=5.17, or on a statin

Co_morbid <- Co_morbid %>%
  mutate(hyperlipidemia_pmh = case_when(
    raised_lipids == 1 | Cholesterol_meds == 1 ~ 1,                       # If either column is 1, set 1
    is.na(raised_lipids) & is.na(Cholesterol_meds) ~ NA_real_,            # If both are NA, set NA
    is.na(raised_lipids) ~ Cholesterol_meds,                              # If only raised_lipids is NA, take cholesterol_meds
    is.na(Cholesterol_meds) ~ raised_lipids,                              # If only cholesterol_meds is NA, take raised_lipids
    TRUE ~ 0                                                              # If both are 0, set 0
  ))

##Hypertension

check_hypertension <- function(text) {
  if (is.na(text)) {
    return(NA)
  } else if (grepl("High blood pressure", text, ignore.case = T)) {
    return(1)
  } else {
    return(0)
  }
}

check_hypertension_self <- function(text) {
  if (is.na(text)) {
    return(NA)
  } else if (grepl("hypertension|essential hypertension", text, ignore.case = T)) {
    return(1)
  } else {
    return(0)
  }
}

Demographics$hypertension_verbal <- sapply(Demographics$p6150_i0, check_hypertension)
Demographics$hypertension_self <- sapply(Demographics$p20002_i0, check_hypertension_self)

Demographics <- Demographics %>%
  mutate(hypertension_recruitment = case_when(
    hypertension_verbal == 1 | hypertension_self == 1 ~ 1,                       # If either column is 1, set 1
    is.na(hypertension_verbal) & is.na(hypertension_self) ~ NA_real_,            # If both are NA, set NA
    is.na(hypertension_verbal) ~ hypertension_self,                              
    is.na(hypertension_self) ~ hypertension_verbal,                              
    TRUE ~ 0                                                              
  ))

Co_morbid$min_systolic <- pmin(Co_morbid$p4080_i0_a0, Co_morbid$p4080_i0_a1, na.rm = TRUE)
Co_morbid$min_diastolic <- pmin(Co_morbid$p4079_i0_a0, Co_morbid$p4079_i0_a1, na.rm = TRUE)
Co_morbid$hypertension_BP_readings <- ifelse(Co_morbid$min_systolic>=140, 1,
                                             ifelse(Co_morbid$min_diastolic>=90, 1, 0 ))

Co_morbid$Hypertension_Hx <- ifelse(!is.na(Co_morbid$p131286) | !is.na(Co_morbid$p131288) | !is.na(Co_morbid$p131290) | 
                                      !is.na(Co_morbid$p131292) | !is.na(Co_morbid$p131294), 1, 0)


Co_morbid <- Co_morbid %>% 
  mutate(earliest_date_hypertension = pmin(p131286,p131288,p131290,p131292,p131294, na.rm = TRUE))

Co_morbid$Time_to_Hypertension <-difftime(Co_morbid$earliest_date_hypertension,Co_morbid$Date_of_Birth_imputed,units=c("weeks"))
Co_morbid$Time_to_Hypertension <-Co_morbid$Time_to_Hypertension/52.1429
Co_morbid$Time_to_Hypertension<-gsub("weeks","",as.character(Co_morbid$Time_to_Hypertension))
Co_morbid$Time_to_Hypertension<-as.numeric(Co_morbid$Time_to_Hypertension)

Co_morbid$Hypertension_Hx_recruit <- ifelse(Co_morbid$Time_to_Hypertension<=Co_morbid$Age_at_recruitment, 1,0)
Co_morbid$Hypertension_Hx_recruit[is.na(Co_morbid$Hypertension_Hx_recruit)] <- 0

Co_morbid <- Co_morbid %>%
  mutate(hypertension_co_morbid_recruit = case_when(
    # If any of the columns is 1, set 1
    Hypertension_Hx_recruit == 1 | hypertension_BP_readings == 1 | Hypertension_meds == 1 ~ 1,
    # If all three are NA, set NA
    is.na(Hypertension_Hx_recruit) & is.na(hypertension_BP_readings) & is.na(Hypertension_meds) ~ NA_real_,
    # If 0 is the highest non-missing value, set 0
    TRUE ~ 0
  ))


##Diabetes

check_diabetes <- function(text) {
  if (is.na(text)) {
    return(NA)
  } else if (grepl("Yes", text, ignore.case = T)) {
    return(1)
  } else {
    return(0)
  }
}

check_diabetes_self <- function(text) {
  if (is.na(text)) {
    return(NA)
  } else if (grepl("diabetes|type 1 diabetes|type 2 diabetes", text, ignore.case = T)) {
    return(1)
  } else {
    return(0)
  }
}

Demographics$diabetes_verbal <- sapply(Demographics$p2443_i0, check_diabetes)
Demographics$diabetes_self <- sapply(Demographics$p20002_i0, check_diabetes_self)

Demographics <- Demographics %>%
  mutate(diabetes_recruitment = case_when(
    diabetes_verbal == 1 | diabetes_self == 1 ~ 1,                       # If either column is 1, set 1
    is.na(diabetes_verbal) & is.na(diabetes_self) ~ NA_real_,            # If both are NA, set NA
    is.na(diabetes_verbal) ~ diabetes_self,                              
    is.na(diabetes_self) ~ diabetes_verbal,                              
    TRUE ~ 0                                                              
  ))


Co_morbid$Diabetes_Hx <- ifelse(!is.na(Co_morbid$p130706) | !is.na(Co_morbid$p130708) | !is.na(Co_morbid$p130710) | 
                                  !is.na(Co_morbid$p130712) | !is.na(Co_morbid$p130714), 1, 0)


Co_morbid <- Co_morbid %>% 
  mutate(earliest_date_diabetes = pmin(p130706,p130708,p130710,p130712,p130714, na.rm = TRUE))

Co_morbid$Time_to_Diabetes <-difftime(Co_morbid$earliest_date_diabetes,Co_morbid$Date_of_Birth_imputed,units=c("weeks"))
Co_morbid$Time_to_Diabetes <-Co_morbid$Time_to_Diabetes/52.1429
Co_morbid$Time_to_Diabetes<-gsub("weeks","",as.character(Co_morbid$Time_to_Diabetes))
Co_morbid$Time_to_Diabetes<-as.numeric(Co_morbid$Time_to_Diabetes)

Co_morbid$Diabetes_co_morbid_recruit <- ifelse(Co_morbid$Time_to_Diabetes<=Co_morbid$Age_at_recruitment, 1,0)

Co_morbid$Diabetes_co_morbid_recruit <- ifelse(Co_morbid$Insulin_meds==1, Co_morbid$Insulin_meds, Co_morbid$Diabetes_co_morbid_recruit)
Co_morbid$Diabetes_co_morbid_recruit[is.na(Co_morbid$Diabetes_co_morbid_recruit)] <- 0

##Cerebrovascular

check_stroke <- function(text) {
  if (is.na(text)) {
    return(NA)
  } else if (grepl("Stroke", text, ignore.case = T)) {
    return(1)
  } else {
    return(0)
  }
}


check_stroke_self <- function(text) {
  if (is.na(text)) {
    return(NA)
  } else if (grepl("subarachnoid haemorrhage|brain haemorrhage|ischaemic stroke", 
                   text, ignore.case = TRUE)) {
    return(1)
  } else {
    return(0)
  }
}

Demographics$stroke_verbal <- sapply(Demographics$p6150_i0, check_stroke)
Demographics$stroke_self <- sapply(Demographics$p20002_i0, check_stroke_self)

Demographics <- Demographics %>%
  mutate(stroke_recruitment = case_when(
    stroke_verbal == 1 | stroke_self == 1 ~ 1,                       # If either column is 1, set 1
    is.na(stroke_verbal) & is.na(stroke_self) ~ NA_real_,            # If both are NA, set NA
    is.na(stroke_verbal) ~ stroke_self,                              
    is.na(stroke_self) ~ stroke_verbal,                              
    TRUE ~ 0                                                              
  ))

Co_morbid$Stroke_Hx <- ifelse(!is.na(Co_morbid$p42006) | !is.na(Co_morbid$p131360) | !is.na(Co_morbid$p131362) | 
                                !is.na(Co_morbid$p131364) | !is.na(Co_morbid$p131366)
                              | !is.na(Co_morbid$p131368)| !is.na(Co_morbid$p131374)| !is.na(Co_morbid$p131376)| !is.na(Co_morbid$p131378), 1, 0)


Co_morbid <- Co_morbid %>% 
  mutate(earliest_date_stroke = pmin(p42006, p131360,p131362,p131364,p131366, p131368,p131374,p131376, p131378, na.rm = TRUE))

Co_morbid$Time_to_Stroke <-difftime(Co_morbid$earliest_date_stroke,Co_morbid$Date_of_Birth_imputed,units=c("weeks"))
Co_morbid$Time_to_Stroke <-Co_morbid$Time_to_Stroke/52.1429
Co_morbid$Time_to_Stroke<-gsub("weeks","",as.character(Co_morbid$Time_to_Stroke))
Co_morbid$Time_to_Stroke<-as.numeric(Co_morbid$Time_to_Stroke)

Co_morbid$Stroke_co_morbid_recruit <- ifelse(Co_morbid$Time_to_Stroke<=Co_morbid$Age_at_recruitment, 1,0)
Co_morbid$Stroke_co_morbid_recruit[is.na(Co_morbid$Stroke_co_morbid_recruit)] <- 0

##Cardiovascular

check_MI <- function(text) {
  if (is.na(text)) {
    return(NA)
  } else if (grepl("Heart attack", text, ignore.case = T)) {
    return(1)
  } else {
    return(0)
  }
}

check_angina <- function(text) {
  if (is.na(text)) {
    return(NA)
  } else if (grepl("Angina", text, ignore.case = T)) {
    return(1)
  } else {
    return(0)
  }
}

check_MI_self <- function(text) {
  if (is.na(text)) {
    return(NA)
  } else if (grepl("heart attack/myocardial infarction", text, ignore.case = T)) {
    return(1)
  } else {
    return(0)
  }
}

check_Angina_self <- function(text) {
  if (is.na(text)) {
    return(NA)
  } else if (grepl("angina", text, ignore.case = T)) {
    return(1)
  } else {
    return(0)
  }
}

Demographics$MI_verbal <- sapply(Demographics$p6150_i0, check_MI)
Demographics$Angina_verbal <- sapply(Demographics$p6150_i0, check_angina)

Demographics$MI_self <- sapply(Demographics$p20002_i0, check_MI_self)
Demographics$Angina_self <- sapply(Demographics$p20002_i0, check_Angina_self)

Demographics <- Demographics %>%
  mutate(Cardiovascular_recruitment = case_when(
    # If any of the columns is 1, set 1
    MI_verbal == 1 | Angina_verbal == 1 | MI_self == 1 | Angina_self == 1 ~ 1,
    # If all four are NA, set NA
    (is.na(Angina_self) & is.na(Angina_verbal)) | (is.na(MI_self) & is.na(MI_verbal)) ~ NA_real_,
    # If 0 is the highest non-missing value, set 0
    TRUE ~ 0
  ))


Co_morbid$Cardiovascular_Hx <- ifelse(!is.na(Co_morbid$p42000) | !is.na(Co_morbid$p131296) 
                                      | !is.na(Co_morbid$p131298) | 
                                        !is.na(Co_morbid$p131300) | !is.na(Co_morbid$p131302)
                                      | !is.na(Co_morbid$p131304)| !is.na(Co_morbid$p131306), 1, 0)


Co_morbid <- Co_morbid %>% 
  mutate(earliest_date_cardiovascular = pmin(p42000, p131296,p131298,p131300,p131302,p131304,p131306, na.rm = TRUE))

Co_morbid$Time_to_Cardiovascular <-difftime(Co_morbid$earliest_date_cardiovascular,Co_morbid$Date_of_Birth_imputed,units=c("weeks"))
Co_morbid$Time_to_Cardiovascular <-Co_morbid$Time_to_Cardiovascular/52.1429
Co_morbid$Time_to_Cardiovascular<-gsub("weeks","",as.character(Co_morbid$Time_to_Cardiovascular))
Co_morbid$Time_to_Cardiovascular<-as.numeric(Co_morbid$Time_to_Cardiovascular)

Co_morbid$Cardiovascular_co_morbid_recruit <- ifelse(Co_morbid$Time_to_Cardiovascular<=Co_morbid$Age_at_recruitment, 1,0)
Co_morbid$Cardiovascular_co_morbid_recruit[is.na(Co_morbid$Cardiovascular_co_morbid_recruit)] <- 0

#Depression

check_depression_self <- function(text) {
  if (is.na(text)) {
    return(NA)
  } else if (grepl("depression|mania/bipolar disorder/manic depression|post-natal depression", 
                   text, ignore.case = TRUE)) {
    return(1)
  } else {
    return(0)
  }
}

Demographics$depression_recruitment <- sapply(Demographics$p20002_i0, check_depression_self)

Co_morbid$Depression_Hx <- ifelse(!is.na(Co_morbid$p130894) | !is.na(Co_morbid$p130896) 
                                  | !is.na(Co_morbid$p130898) | 
                                    !is.na(Co_morbid$p130900) | !is.na(Co_morbid$p130902), 1, 0)


Co_morbid <- Co_morbid %>% 
  mutate(earliest_date_depression = pmin(p130894, p130896,p130898, p130900, p130902, na.rm = TRUE))

Co_morbid$Time_to_Depression <-difftime(Co_morbid$earliest_date_depression,Co_morbid$Date_of_Birth_imputed,units=c("weeks"))
Co_morbid$Time_to_Depression <-Co_morbid$Time_to_Depression/52.1429
Co_morbid$Time_to_Depression<-gsub("weeks","",as.character(Co_morbid$Time_to_Depression))
Co_morbid$Time_to_Depression<-as.numeric(Co_morbid$Time_to_Depression)

Co_morbid$Depression_co_morbid_recruit <- ifelse(Co_morbid$Time_to_Depression<=Co_morbid$Age_at_recruitment, 1,0)
Co_morbid$Depression_co_morbid_recruit[is.na(Co_morbid$Depression_co_morbid_recruit)] <- 0

Comorbid_refined <- Co_morbid %>% select(eid, Depression_co_morbid_recruit,Cardiovascular_co_morbid_recruit,Stroke_co_morbid_recruit,Diabetes_co_morbid_recruit,hypertension_co_morbid_recruit, hyperlipidemia_pmh)

Demographics <- merge(Demographics, Comorbid_refined, by="eid", all.x = T)

Demographics <- Demographics %>%
  mutate(Hypertension_pmh = case_when(
    hypertension_co_morbid_recruit == 1 | hypertension_recruitment == 1 ~ 1,                       
    is.na(hypertension_co_morbid_recruit) & is.na(hypertension_recruitment) ~ NA_real_,            
    is.na(hypertension_co_morbid_recruit) ~ hypertension_recruitment,                              
    is.na(hypertension_recruitment) ~ hypertension_co_morbid_recruit,                              
    TRUE ~ 0                                                             
  ))

Demographics <- Demographics %>%
  mutate(Diabetes_pmh = case_when(
    Diabetes_co_morbid_recruit == 1 | diabetes_recruitment == 1 ~ 1,                       
    is.na(Diabetes_co_morbid_recruit) & is.na(diabetes_recruitment) ~ NA_real_,            
    is.na(Diabetes_co_morbid_recruit) ~ diabetes_recruitment,                              
    is.na(diabetes_recruitment) ~ Diabetes_co_morbid_recruit,                              
    TRUE ~ 0                                                              
  ))


Demographics <- Demographics %>%
  mutate(Cerebrovascular_pmh = case_when(
    Stroke_co_morbid_recruit == 1 | stroke_recruitment == 1 ~ 1,                       
    is.na(Stroke_co_morbid_recruit) & is.na(stroke_recruitment) ~ NA_real_,            
    is.na(Stroke_co_morbid_recruit) ~ stroke_recruitment,                              
    is.na(stroke_recruitment) ~ Stroke_co_morbid_recruit,                             
    TRUE ~ 0                                                             
  ))

Demographics <- Demographics %>%
  mutate(Cardiovascular_pmh = case_when(
    Cardiovascular_co_morbid_recruit == 1 | Cardiovascular_recruitment == 1 ~ 1,                      
    is.na(Cardiovascular_co_morbid_recruit) & is.na(Cardiovascular_recruitment) ~ NA_real_,           
    is.na(Cardiovascular_co_morbid_recruit) ~ Cardiovascular_recruitment,                              
    is.na(Cardiovascular_recruitment) ~ Cardiovascular_co_morbid_recruit,                              
    TRUE ~ 0                                                              
  ))

Demographics <- Demographics %>%
  mutate(Depression_pmh = case_when(
    Depression_co_morbid_recruit == 1 | depression_recruitment == 1 ~ 1,                      
    is.na(Depression_co_morbid_recruit) & is.na(depression_recruitment) ~ NA_real_,            
    is.na(Depression_co_morbid_recruit) ~ depression_recruitment,                             
    is.na(depression_recruitment) ~ Depression_co_morbid_recruit,                              
    TRUE ~ 0                                                              
  ))

###########Cleaning the head injured cohort##################

Head_injured$age_category_1 <- ifelse(is.na(Head_injured$age_category_1), "unknown", Head_injured$age_category_1)

##Finding those with only GP recording of first or second injury

Head_injured$GP_only <- ifelse(!Head_injured$source_1=="Community",0,1)
Head_injured$HES_source <- ifelse(!Head_injured$source_1=="HES",0,1)

#26925 head injured

age_at_death <- Death_reg %>% select(eid, p40007_i0)
Head_injured <- merge(Head_injured, age_at_death, by = "eid", all.x = T)
Head_injured$Time_to_injury_1 <- as.numeric(Head_injured$Time_to_injury_1)
Head_injured$p40007_i0 <- as.numeric(Head_injured$p40007_i0)

Head_injured_at_death <- subset(Head_injured, abs(Time_to_injury_1 - p40007_i0) <= 0.25)


Head_injured <- subset(Head_injured, !Head_injured$eid %in% Head_injured_at_death$eid)


Head_injured <- Head_injured %>% select(-c(p40007_i0))
#26240 head injured

Head_injured_refined <- Head_injured %>% select(c("eid", "source_1", "source_2", "Time_to_injury_1", "Time_to_injury_2", "GP_only", "age_category_1", "age_category_2", "HES_source", "severity_1"))
Head_injured_refined$cohort <- "Yes"
Head_injured$cohort <- "Yes"

#Overall <- merge(Head_injured_refined, Demographics, by="eid", all.y = T)

Overall <- merge(Head_injured, Demographics, by="eid", all.y = T)
Nero <- Nero %>% select(-c(Date_of_Birth_imputed))
Overall <- merge(Nero, Overall, by="eid", all.y = T)

#Overall$date_of_HES_censoring <- ifelse(
#Overall$p40022 == "Hospital Episode Statistics from England", "2023-11-30",
# ifelse(Overall$p40022 == "Patient Episode Database for Wales", "2023-11-30",
#     ifelse(Overall$p40022 == "Scottish Morbidity Records", "2023-12-31", NA)))

#Overall$date_of_HES_censoring[is.na(Overall$p40022)] <- "2023-11-30"

Overall$date_of_HES_censoring <- "2022-10-31"

Overall$date_of_HES_censoring<-as.factor(Overall$date_of_HES_censoring)
Overall$date_of_HES_censoring<-strptime(Overall$date_of_HES_censoring,format="%Y-%m-%d")
Overall$date_of_HES_censoring<-as.Date(Overall$date_of_HES_censoring,format="%Y-%m-%d")

Overall$earliest_date_neurodegeneration<-as.factor(Overall$earliest_date_neurodegeneration)
Overall$earliest_date_neurodegeneration<-strptime(Overall$earliest_date_neurodegeneration,format="%Y-%m-%d")
Overall$earliest_date_neurodegeneration<-as.Date(Overall$earliest_date_neurodegeneration,format="%Y-%m-%d")

##Ignoring those with the outcome after the date of censoring

Neurodegen_after_cut <- subset(Overall, Overall$date_of_HES_censoring < Overall$earliest_date_neurodegeneration)
pre <- subset(Overall, Overall$neurodegeneration==1)
check3 <- Overall
columns_to_modify <- c(
  "all_cause_dementia", "Alzhiemers", "Vascular", "Frontotemporal", 
  "Time_to_all_cause_dementia", "Time_to_Alzhiemers", 
  "Time_to_Vasc", "Time_to_FTD", "Parkinsons_present",  
  "Time_to_PD", "MND_present",  "Time_to_MND", 
  "neurodegeneration",  "Time_to_Neurodegeneration"
)

Overall[columns_to_modify] <- lapply(Overall[columns_to_modify], function(column) {
  ifelse(Overall$eid %in% Neurodegen_after_cut$eid, NA, column)
})

check2 <- subset(Overall, Overall$eid %in% Neurodegen_after_cut$eid)
post <- subset(Overall, Overall$neurodegeneration==1)
Overall$p191<-as.factor(Overall$p191)
Overall$p191<-strptime(Overall$p191,format="%Y-%m-%d")
Overall$p191<-as.Date(Overall$p191,format="%Y-%m-%d")

Overall$Date_of_Birth_imputed<-as.factor(Overall$Date_of_Birth_imputed)
Overall$Date_of_Birth_imputed<-strptime(Overall$Date_of_Birth_imputed,format="%Y-%m-%d")
Overall$Date_of_Birth_imputed<-as.Date(Overall$Date_of_Birth_imputed,format="%Y-%m-%d")

Overall$Age_at_HES_censoring <- difftime(Overall$date_of_HES_censoring,Overall$Date_of_Birth_imputed,units=c("weeks"))
Overall$Age_at_HES_censoring <-Overall$Age_at_HES_censoring/52.1429
Overall$Age_at_HES_censoring<-gsub("weeks","",as.character(Overall$Age_at_HES_censoring))
Overall$Age_at_HES_censoring<-as.numeric(Overall$Age_at_HES_censoring)

##Ignoring those with head injury after the date of censoring

Head_inj_after_cut <- subset(Overall, Overall$Age_at_HES_censoring < Overall$Time_to_injury_1)
pre <- subset(Overall, Overall$cohort=="Yes")
check3 <- Overall
columns_to_modify <- c(
  "source_1", "source_2", "Time_to_injury_1", "Time_to_injury_2", "GP_only", "age_category_1", "age_category_2", "severity_1"
)

Overall[columns_to_modify] <- lapply(Overall[columns_to_modify], function(column) {
  ifelse(Overall$eid %in% Head_inj_after_cut$eid, NA, column)
})

check2 <- subset(Overall, Overall$eid %in% Head_inj_after_cut$eid)
post <- subset(Overall, Overall$cohort=="Yes")


modify_columns <- function(Overall, suffixes, cohort_condition = "Yes") {
  Head_inj_after_cut <- subset(Overall, Overall$Age_at_HES_censoring < Overall[[paste0("Time_to_injury_", suffixes[1])]])
  pre <- subset(Overall, Overall$cohort == cohort_condition)
  columns_to_modify <- paste0(c("source", "Time_to_injury", "age_category", "severity"), "_", suffixes)
  Overall[columns_to_modify] <- lapply(Overall[columns_to_modify], function(column) {
    ifelse(Overall$eid %in% Head_inj_after_cut$eid, NA, column)
  })
  
  return(Overall)
}

suffixes <- c(2, 3,4,5,6,7,8)  # Can modify this list with any suffixes you want
modified_Overall <- modify_columns(Overall, suffixes)

Head_inj_after_cut <- subset(Overall, Overall$Age_at_HES_censoring < Overall$Time_to_injury_2)
pre <- subset(Overall, Overall$cohort=="Yes")
check3 <- Overall
columns_to_modify <- c(
  "source_2",  "Time_to_injury_2", "age_category_2"
)

Overall[columns_to_modify] <- lapply(Overall[columns_to_modify], function(column) {
  ifelse(Overall$eid %in% Head_inj_after_cut$eid, NA, column)
})

check2 <- subset(Overall, Overall$eid %in% Head_inj_after_cut$eid)
post <- subset(Overall, Overall$cohort=="Yes")


Overall$Age_at_loss_to_follow <- difftime(Overall$p191,Overall$Date_of_Birth_imputed,units=c("weeks"))
Overall$Age_at_loss_to_follow <-Overall$Age_at_loss_to_follow/52.1429
Overall$Age_at_loss_to_follow<-gsub("weeks","",as.character(Overall$Age_at_loss_to_follow))
Overall$Age_at_loss_to_follow<-as.numeric(Overall$Age_at_loss_to_follow)



Overall$age_at_death_or_HES <- ifelse(
  is.na(Overall$p40007_i0),
  Overall$Age_at_HES_censoring,
  ifelse(Overall$p40007_i0<Overall$Age_at_HES_censoring,
         Overall$p40007_i0,
         Overall$Age_at_HES_censoring)
)

Overall$Age_at_loss_to_follow <- ifelse(
  !is.na(Overall$Time_to_injury_1) & 
    !is.na(Overall$Age_at_loss_to_follow) & 
    Overall$Time_to_injury_1 > Overall$Age_at_loss_to_follow,
  Overall$age_at_death_or_HES,
  Overall$Age_at_loss_to_follow
)

Overall$Age_at_loss_to_follow <- ifelse(
  !is.na(Overall$Time_to_Neurodegeneration) & 
    !is.na(Overall$Age_at_loss_to_follow) & 
    Overall$Time_to_Neurodegeneration > Overall$Age_at_loss_to_follow,
  Overall$age_at_death_or_HES,
  Overall$Age_at_loss_to_follow
)

Overall$Age_at_censoring <- ifelse(
  is.na(Overall$Age_at_loss_to_follow),
  Overall$age_at_death_or_HES,
  ifelse(Overall$Age_at_loss_to_follow<Overall$age_at_death_or_HES,
         Overall$Age_at_loss_to_follow,
         Overall$age_at_death_or_HES)
)

Overall <- Overall %>%
  select(-starts_with("X"))

Overall <- Overall %>%
  mutate(
    greatest_severity = case_when(
      is.na(severity_1) ~ NA_character_,  # Set to NA if severity_1 is NA
      rowSums(select(., starts_with("severity_")) == "moderate/severe", na.rm = TRUE) > 0 ~ "moderate/severe",
      TRUE ~ "mild"
    )
  )

check <- subset(Overall, Overall$Time_to_injury_1>Overall$Age_at_censoring)
check <- check %>% select(c(eid, Time_to_injury_1, Age_at_censoring, p190, p191, Date_of_Birth_imputed, Age_at_loss_to_follow, age_at_death_or_HES, p40007_i0))

####################Creating the analysis dataframe##############################################

For_imputation <- Overall
For_imputation <- For_imputation %>%
  mutate_all(~replace(., . %in% c("Do not know", "Prefer not to answer"), NA))

##Step 1: Remove any subjects with neurodegeneration prior to recruitment to the Biobank, or if date/time is NA

date_columns <- c("p53_i0")

for (col in date_columns) {
  
  Overall[[col]] <- as.factor(Overall[[col]])
  
  Overall[[col]] <- strptime(Overall[[col]], format="%Y-%m-%d")
  
  Overall[[col]] <- as.Date(Overall[[col]], format="%Y-%m-%d")
}

Overall$prior_neurodegen <- as.integer(Overall$p53_i0 >= Overall$earliest_date_neurodegeneration)

Overall$prior_neurodegen <- ifelse((Overall$neurodegeneration==1&is.na(Overall$earliest_date_neurodegeneration)),1,Overall$prior_neurodegen)

Overall <- subset(Overall, Overall$prior_neurodegen==0 | is.na(Overall$prior_neurodegen))

##Step 2: Remove Neonatal, paediatric head injury or people with no time to head injury

table(Overall$age_category_1) #500480
Overall <- subset(Overall, is.na(Overall$source_NA)) #500465
Overall <- subset(Overall, Overall$age_category_1=="Adult" | is.na(Overall$age_category_1)) #498779

##Step 3: Exclude with no birth data

Overall <- subset(Overall, !is.na(Overall$Date_of_Birth_imputed))

For_biom <- Overall
For_biom <- For_biom %>%
  mutate_all(~replace(., . %in% c("Do not know", "Prefer not to answer"), NA))

####################Demographics table############################################################

##Do not know or prefer not to answer to NA
Overall <- Overall %>%
  mutate_all(~replace(., . %in% c("Do not know", "Prefer not to answer"), NA))
Overall$p21022 <- as.numeric(Overall$p21022)
Overall$p22189 <- as.numeric(Overall$p22189)
Overall$p21001_i0 <- as.numeric(Overall$p21001_i0)

Overall$highest_level <- factor(Overall$highest_level, levels = c("College", "Training", "Secondary", "Primary"))
Overall$p2178 <- factor(Overall$p2178, levels = c("Excellent", "Good", "Fair", "Poor", "Prefer not to answer", "Do not know"))
Overall$p2178 <- droplevels(Overall$p2178)
Overall$p21000_i0 <- factor(Overall$p21000_i0, levels = c("Asian British", "Black British", "White British", "Other"))
Overall$p21000_i0 <- droplevels(Overall$p21000_i0)
Overall$p2020_i0 <- factor(
  Overall$p2020_i0,
  levels = c( "No", "Yes", "Prefer not to answer", "Do not know")
)
Overall$p2020_i0 <- droplevels(Overall$p2020_i0)

Overall$p22032_i0 <- factor(
  Overall$p22032_i0,
  levels = c( "high", "moderate", "low")
)
Overall$p22032_i0 <- droplevels(Overall$p22032_i0)
Overall$Lives_alone <- factor(
  Overall$Lives_alone,
  levels = c("Lives Alone", "Lives with others", "Prefer not to answer", "Do not know")
)
Overall$Lives_alone <- droplevels(Overall$Lives_alone)

Head_injury_only <- subset(Overall, Overall$cohort=="Yes")

Overall$cohort <- ifelse(is.na(Overall$cohort),"No", Overall$cohort)


####What type of neurodegeneration

Overall$all_cause_dementia <- ifelse(is.na(Overall$all_cause_dementia), 0, Overall$all_cause_dementia)
Overall$Alzhiemers <- ifelse(is.na(Overall$Alzhiemers), 0, Overall$Alzhiemers)
Overall$Vascular <- ifelse(is.na(Overall$Vascular), 0, Overall$Vascular)
Overall$Frontotemporal <- ifelse(is.na(Overall$Frontotemporal), 0, Overall$Frontotemporal)
Overall$Parkinsons <- ifelse(is.na(Overall$Parkinsons_present), 0, Overall$Parkinsons_present)
Overall$MND <- ifelse(is.na(Overall$MND_present), 0, Overall$MND_present)
Overall$neurodegeneration <- ifelse(is.na(Overall$neurodegeneration), 0, Overall$neurodegeneration)


###################Main analysis###################################################

##First simplify the Outcome_dataset##

Overall$Neurogeneration_diagnosis <- Overall$neurodegeneration

Overall$Neurogeneration_diagnosis[is.na(Overall$Neurogeneration_diagnosis)] <- 0

Overall$dead <- ifelse(!is.na(Overall$p40007_i0), 1, 0)
Outcome_dataset <- Overall %>% select(c(eid, Date_of_Birth_imputed,dead, Age_at_censoring,p40007_i0, p53_i0, cohort, Neurogeneration_diagnosis,p40007_i0, Age_at_HES_censoring, Time_to_Neurodegeneration, Time_to_injury_1, Time_to_injury_2, Age_at_recruitment, Time_to_all_cause_dementia, all_cause_dementia, Time_to_PD, Parkinsons, GP_only, MND))

#To create fine gray model outcome for competing risks you need a variable that is 0 if censored, 1 if neurodegen and 2 if dead. Therefore I can code for a 1 if 
#1 is in the neurodegen diagnosis, then indicate people who died after neurodegen (these will be NA)

Outcome_dataset$dead_after_neurodegeneration <- ifelse(
  !is.na(Outcome_dataset$p40007_i0) & Outcome_dataset$Time_to_Neurodegeneration <= Outcome_dataset$p40007_i0, 
  1, 
  NA
)

Outcome_dataset$Fine_gray_outcome <- with(Outcome_dataset, ifelse(
  Neurogeneration_diagnosis == 1, 1,
  ifelse(dead == 1 & (is.na(p40007_i0) | is.na(dead_after_neurodegeneration)), 2, 0)
))

#As the date of birth was imputed, some died re recruitment, change this

Outcome_dataset$Age_at_censoring <- ifelse(Outcome_dataset$Age_at_censoring < Outcome_dataset$Age_at_recruitment, (Outcome_dataset$Age_at_recruitment+0.01), Outcome_dataset$Age_at_censoring)

# Calculate the end time considering neurodegeneration and censoring
Outcome_dataset$end <- ifelse(is.na(Outcome_dataset$Time_to_Neurodegeneration), 
                              Outcome_dataset$Age_at_censoring, 
                              Outcome_dataset$Time_to_Neurodegeneration)
Outcome_dataset$end <- as.numeric(Outcome_dataset$end)


Outcome_dataset$head_injury_start_1 <- ifelse(is.na(Outcome_dataset$Time_to_injury_1), NA, Outcome_dataset$Time_to_injury_1)
Outcome_dataset$head_injury_start_1 <- ifelse(Outcome_dataset$head_injury_start_1 < Outcome_dataset$Age_at_recruitment, 
                                              (Outcome_dataset$Age_at_recruitment+0.01), 
                                              Outcome_dataset$head_injury_start_1)

Outcome_dataset$head_injury_start_1 <- Outcome_dataset$head_injury_start_1 - Outcome_dataset$Age_at_recruitment

Outcome_dataset$head_injury_start_2 <- ifelse(is.na(Outcome_dataset$Time_to_injury_2), NA, Outcome_dataset$Time_to_injury_2)
Outcome_dataset$head_injury_start_2 <- ifelse(Outcome_dataset$head_injury_start_2 < Outcome_dataset$Age_at_recruitment, 
                                              (Outcome_dataset$Age_at_recruitment+0.015), 
                                              Outcome_dataset$head_injury_start_2)

Outcome_dataset$head_injury_start_2 <- Outcome_dataset$head_injury_start_2 - Outcome_dataset$Age_at_recruitment

Outcome_dataset$start <- 0
Outcome_dataset$stop <- Outcome_dataset$end - Outcome_dataset$Age_at_recruitment
Outcome_dataset$event <- Outcome_dataset$Neurogeneration_diagnosis
Outcome_dataset$time_in_study <- Outcome_dataset$stop

#####Dementia###

Outcome_dataset$all_cause_dementia[is.na(Outcome_dataset$all_cause_dementia)] <- 0

Outcome_dataset$end_dementia <- ifelse(is.na(Outcome_dataset$Time_to_all_cause_dementia), 
                                       Outcome_dataset$Age_at_censoring, 
                                       Outcome_dataset$Time_to_all_cause_dementia)
Outcome_dataset$end_dementia <- as.numeric(Outcome_dataset$end_dementia)

Outcome_dataset$stop_dementia <- Outcome_dataset$end_dementia - Outcome_dataset$Age_at_recruitment
Outcome_dataset$event_dementia <- Outcome_dataset$all_cause_dementia
Outcome_dataset$time_in_study_dementia <- Outcome_dataset$stop_dementia

#############Parkinsons###########################################

Outcome_dataset$Parkinsons[is.na(Outcome_dataset$Parkinsons)] <- 0


Outcome_dataset$end_Parkinsons <- ifelse(is.na(Outcome_dataset$Time_to_PD), 
                                         Outcome_dataset$Age_at_censoring, 
                                         Outcome_dataset$Time_to_PD)
Outcome_dataset$end_Parkinsons <- as.numeric(Outcome_dataset$end_Parkinsons)

Outcome_dataset$stop_Parkinsons <- Outcome_dataset$end_Parkinsons - Outcome_dataset$Age_at_recruitment
Outcome_dataset$event_Parkinsons <- Outcome_dataset$Parkinsons 
Outcome_dataset$time_in_study_Parkinsons <- Outcome_dataset$stop_Parkinsons

Simple <- Outcome_dataset %>% select(c(eid, start, stop, event, head_injury_start_1, head_injury_start_2, time_in_study, stop_dementia, event_dementia, time_in_study_dementia, stop_Parkinsons, event_Parkinsons, time_in_study_Parkinsons, GP_only, Parkinsons, MND, all_cause_dementia))

##three months lag ########################################################################

Outcome_dataset_3_months_lag <- Outcome_dataset

Outcome_dataset_3_months_lag$Time_to_injury_1 <- as.numeric(Outcome_dataset_3_months_lag$Time_to_injury_1)
Outcome_dataset_3_months_lag$Time_to_Neurodegeneration <- as.numeric(Outcome_dataset_3_months_lag$Time_to_Neurodegeneration)

initial_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_3_months_lag[[time_col]]))
})

initial_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("Time_to_injury_", i, sep = "")
  sum(is.na(Outcome_dataset_3_months_lag[[time_col]]))
})

for (i in 1:2) {
  start_col <- paste("head_injury_start_", i, sep = "")
  gg_col <- paste("Time_to_injury_", i, sep = "")
  time_col <- paste("Time_to_injury_", i, sep = "")
  condition <- Outcome_dataset_3_months_lag[[time_col]] > (Outcome_dataset_3_months_lag$Time_to_Neurodegeneration - 0.25)
  Outcome_dataset_3_months_lag[[start_col]][condition] <- NA
  Outcome_dataset_3_months_lag[[gg_col]][condition] <- NA
}

updated_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_3_months_lag[[time_col]]))
})

updated_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("Time_to_injury_", i, sep = "")
  sum(is.na(Outcome_dataset_3_months_lag[[time_col]]))
})


na_changes <- updated_na_counts - initial_na_counts
print(na_changes)


Simple_3_months_lag <- Outcome_dataset_3_months_lag %>% select(c(eid, start, stop, event, Fine_gray_outcome, head_injury_start_1, head_injury_start_2, time_in_study, stop_dementia, event_dementia, time_in_study_dementia, stop_Parkinsons, event_Parkinsons, time_in_study_Parkinsons, GP_only, Parkinsons, MND, all_cause_dementia, Time_to_injury_1))


########### 1 year lag  ##############################################################################################################

Outcome_dataset_1_year_lag <- Outcome_dataset

Outcome_dataset_1_year_lag$Time_to_injury_1 <- as.numeric(Outcome_dataset_1_year_lag$Time_to_injury_1)
Outcome_dataset_1_year_lag$Time_to_Neurodegeneration <- as.numeric(Outcome_dataset_1_year_lag$Time_to_Neurodegeneration)

initial_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_1_year_lag[[time_col]]))
})

for (i in 1:2) {
  start_col <- paste("head_injury_start_", i, sep = "")
  time_col <- paste("Time_to_injury_", i, sep = "")
  condition <- Outcome_dataset_1_year_lag[[time_col]] > (Outcome_dataset_1_year_lag$Time_to_Neurodegeneration - 1)
  Outcome_dataset_1_year_lag[[start_col]][condition] <- NA
}

updated_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_1_year_lag[[time_col]]))
})

# Compare before and after
na_changes <- updated_na_counts - initial_na_counts
print(na_changes)

Simple_1_year_lag <- Outcome_dataset_1_year_lag %>% select(c(eid, start, stop, event, head_injury_start_1, head_injury_start_2, time_in_study, stop_dementia, event_dementia, time_in_study_dementia, stop_Parkinsons, event_Parkinsons, time_in_study_Parkinsons, GP_only, Parkinsons, MND, all_cause_dementia))

##########3 year lag

Outcome_dataset_3_year_lag <- Outcome_dataset

Outcome_dataset_3_year_lag$Time_to_injury_1 <- as.numeric(Outcome_dataset_3_year_lag$Time_to_injury_1)
Outcome_dataset_3_year_lag$Time_to_Neurodegeneration <- as.numeric(Outcome_dataset_3_year_lag$Time_to_Neurodegeneration)

initial_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_3_year_lag[[time_col]]))
})

for (i in 1:2) {
  start_col <- paste("head_injury_start_", i, sep = "")
  time_col <- paste("Time_to_injury_", i, sep = "")
  condition <- Outcome_dataset_3_year_lag[[time_col]] > (Outcome_dataset_3_year_lag$Time_to_Neurodegeneration - 3)
  Outcome_dataset_3_year_lag[[start_col]][condition] <- NA
}

updated_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_1_year_lag[[time_col]]))
})


na_changes <- updated_na_counts - initial_na_counts
print(na_changes)

Simple_3_year_lag <- Outcome_dataset_3_year_lag %>% select(c(eid, start, stop, event, head_injury_start_1, head_injury_start_2, time_in_study, stop_dementia, event_dementia, time_in_study_dementia, stop_Parkinsons, event_Parkinsons, time_in_study_Parkinsons, GP_only, Parkinsons, MND, all_cause_dementia))


####Cox Proportional Hazards model#############

###Non adjusted

##categorisation


Overall$cohort <- ifelse(!is.na(Overall$Time_to_injury_1),"Yes", "No")
Overall$neuro_pre <- ifelse(Overall$Time_to_Neurodegeneration<Overall$Time_to_injury_1+0.25, 1, 0)
Overall$cohort_of_interest <- ifelse(!is.na(Overall$neuro_pre) & Overall$neuro_pre == 1, "No", Overall$cohort)

Non<- c("p21022","p22189","p21001_i0")
myVars <- c("p21022", "p31", "p22189", "highest_level", "smoking", "alcohol","p21000_i0","p22032_i0",  "p21001_i0","p2178", "Lives_alone", "p2020_i0")
catVars <- c("p31", "highest_level","p20116_i0", "p21000_i0","p22032_i0",  "p2178", "Lives_alone", "p2020_i0", "smoking", "alcohol")

tab <- CreateTableOne(vars = myVars, data = Overall, strata="cohort_of_interest", factorVars = catVars)
tab <- print(tab, formatOptions = list(big.mark = ","), nonnormal =Non, includeNa = T)
tab
write.csv(tab, "Neurodegeneration_table_one_TBI.csv")

Non<- c("Time_to_all_cause_dementia","Time_to_PD","Time_to_MND")
myVars <- c("neurodegeneration", "all_cause_dementia","Alzhiemers","Vascular","Frontotemporal","Time_to_all_cause_dementia", "Parkinsons", "Time_to_PD", "MND", "Time_to_MND")
catVars <- c("neurodegeneration","all_cause_dementia","Alzhiemers","Vascular","Frontotemporal",  "Parkinsons","MND")

tab <- CreateTableOne(vars = myVars, data = Overall,factorVars = catVars)
tab <- print(tab, formatOptions = list(big.mark = ","), nonnormal =Non, includeNa = T)
tab
write.csv(tab, "Neurodegeneration_type_TBI.csv")

tab <- CreateTableOne(vars = myVars, data = Overall,strata="cohort_of_interest",factorVars = catVars)
tab <- print(tab, formatOptions = list(big.mark = ","), nonnormal =Non, includeNa = T)

write.csv(tab, "Neurodegeneration_type_stratified_TBI.csv")

tab <- CreateTableOne(vars = myVars, data = Overall, factorVars = catVars)
tab <- print(tab, formatOptions = list(big.mark = ",", scientific = FALSE), nonnormal =Non, includeNa = T)
tab
write.csv(tab, "Neurodegeneration_table_one_overall_TBI.csv")

####################################3 month lag

td_df <- tmerge(Simple_3_months_lag, Simple_3_months_lag, id = eid, outcome = event(stop, event))
td_df <- tmerge(td_df, Simple_3_months_lag, id = eid, head_injury = tdc(head_injury_start_1))
td_df$head_injury <- as.factor(td_df$head_injury)
cox_fit_A <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = td_df, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_A)

sum(td_df$tstop-td_df$tstart)
sum(Simple_1_year_lag$stop)

export_cox_summary_to_csv(cox_fit_A, "cox_fit_Unadjusted_TBI.csv")

crude_incidence <- td_df %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )


#######Kapplain myer########

s1 <- survfit(Surv(tstart, tstop, outcome) ~  head_injury, data = td_df, id = eid, timefix = TRUE)
summary_s1 <- summary(s1)

x_start <- min(td_df$tstart, na.rm = TRUE)
x_stop <- max(td_df$tstop, na.rm = TRUE)

Kaplain_Myer <- ggsurvplot(s1, fun = "event", 
                           xlim = c(0, 15),  # Adjust the plot x-axis to 15 years
                           ylim = c(0, 0.1), 
                           risk.table = TRUE, 
                           risk.table.xlim = c(0, 15),  # Align the risk table to 15 years
                           conf.int = TRUE, 
                           break.x.by = 5,  # Major breaks every 5 years
                           legend.labs = c("No Head Injury", "Head Injury"),
                           ggtheme = theme_classic() +
                             theme(
                               legend.text = element_text(size = 12),  
                               legend.title = element_blank(), 
                               axis.title.x = element_text(size = 14),  
                               axis.title.y = element_text(size = 14),
                               axis.text = element_text(size=12),
                               plot.title = element_text(size = 18),
                               legend.position = "right"  ),
                           title = "",
                           xlab = "Time since recruitment (years)", 
                           ylab = "Cumulative neurodegeneration" ,
                           censor = FALSE 
                           
                           )
Kaplain_Myer$plot <- Kaplain_Myer$plot + guides(color = guide_legend(ncol = 1))
Kaplain_Myer$plot <- Kaplain_Myer$plot +
  scale_x_continuous(limits = c(0, 15), breaks = c(0, 5, 10, 15))

ggsave("cumulative_incidence_plot_TBI.tif", plot = Kaplain_Myer$plot, width = 10, height = 6)
ggsave("cumulative_incidence_risk_table_TBI.tif", plot = Kaplain_Myer$table, width = 10, height = 6)

###multiple head injuries


td_df_multi <- tmerge(Simple_3_months_lag, Simple_3_months_lag, id = eid, outcome = event(stop, event))


td_df_multi <- tmerge(td_df_multi, Simple_3_months_lag, id = eid, 
                      head_injury = tdc(head_injury_start_1))


td_df_multi <- tmerge(td_df_multi, Simple_3_months_lag, id = eid, 
                      second_head_injury = tdc(head_injury_start_2))

td_df_multi <- td_df_multi %>%
  mutate(
    head_injury_status = case_when(
      head_injury == 0 ~ "No head injury",
      head_injury == 1 & (is.na(second_head_injury) | second_head_injury == 0) ~ "1 head injury",
      head_injury == 1 & second_head_injury == 1 ~ "2+ head injuries"
    ),
    head_injury_status = factor(head_injury_status, levels = c("No head injury", "1 head injury", "2+ head injuries"))
  )

cox_fit <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury_status, 
                 data = td_df_multi, id = eid, ties = 'efron', timefix = TRUE)

export_cox_summary_to_csv(cox_fit, "cox_fit_Unadjusted_multiple_TBI.csv")
summary(cox_fit)

s1 <- survfit(Surv(tstart, tstop, outcome) ~ head_injury_status, data = td_df_multi, id = eid, timefix = TRUE)

# Define the x-axis limits based on the minimum and maximum of tstart and tstop
x_start <- min(td_df_multi$tstart, na.rm = TRUE)
x_stop <- max(td_df_multi$tstop, na.rm = TRUE)

Kaplan_Meyer_plot_multi <- ggsurvplot(
  s1, 
  fun = "event", 
  xlim = c(0, 15),  # Adjust the plot x-axis to 15 years
  ylim = c(0, 0.2), 
  risk.table = TRUE, 
  risk.table.xlim = c(0, 15),  
  conf.int = TRUE, 
  break.x.by = 5,  # Major breaks every 5 years
  legend.labs = c("No Head Injury", "1 Head Injury", "2+ Head Injuries"),
  ggtheme = theme_classic() +
    theme(
      legend.text = element_text(size = 12),  
      legend.title = element_blank(), 
      axis.title.x = element_text(size = 14),  
      axis.title.y = element_text(size = 14),
      axis.text = element_text(size=12),
      plot.title = element_text(size = 18),
      legend.position = "right"  ),     
  title = "",
  xlab = "Time since recruitment (years)", 
  ylab = "Cumulative neurodegeneration",
  censor = FALSE
)
Kaplan_Meyer_plot_multi$plot <- Kaplan_Meyer_plot_multi$plot + guides(color = guide_legend(ncol = 1))
Kaplan_Meyer_plot_multi$plot <- Kaplan_Meyer_plot_multi$plot +
  scale_x_continuous(limits = c(0, 15), breaks = c(0, 5, 10, 15))

# Display the updated Kaplan-Meier plot
Kaplan_Meyer_plot_multi


# Save the cumulative incidence plot and risk table as images
ggsave("cumulative_incidence_plot_multiple_TBI.tif", plot = Kaplan_Meyer_plot_multi$plot, width = 10, height = 6)
ggsave("cumulative_incidence_risk_table_multiple_TBI.tif", plot = Kaplan_Meyer_plot_multi$table, width = 10, height = 6)

combined_plot <- plot_grid(
  Kaplain_Myer$plot,  
  Kaplan_Meyer_plot_multi$plot,  
  labels = "AUTO",
  label_size = 30,
  ncol = 1,  # Arrange vertically
  rel_heights = c(1, 1)  # Bottom plot is 1.5 times bigger than the top plot
)

# Save the combined plot
ggsave("combined_Kaplan_Meyer_plot.tif", plot = combined_plot, width = 10, height = 12)

# Save the combined image
ggsave(
  filename = "Kaplan_Meyer_plot_with_table_TBI.tif",
  plot = combined_plot,
  width = 10,    # Adjust the width
  height = 6,   # Adjust the height to accommodate the table
  dpi = 300      # High resolution
)



#############Multivariable    

#Age at recruitment, Alcohol (alcohol), Depression (Depression_pmh), Deprivatio 22189), Diabetes (Diabetes_pmh), Education Level (highest_level), Family history (FHDem, FHPD), Genetic Predisposition (p26206, p26260), Hyperlipidemia (hyperlipidemia_pmh), Hypertension (Hypertension_pmh), Obesity (obesity), 
#Physical Activity (p22032_i0), Sex (p31), Smoking (smoking), Social Isolation (p2020_i0), hearing_loss, visual_loss, airpollution(p24003, p24006) (Age_at_recruitment+alcohol+Depression_pmh+Diabetes_pmh+highest_level+DemFH+PDFH+p26206+p26260+hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0)

Adj_1_cov <- Overall %>% select(c(eid, Age_at_recruitment, p31, p2020_i0, p22189, highest_level, p24003, p24006, p22032_i0, obesity, smoking,alcohol))

Missing <- Overall %>% select(c(Age_at_recruitment, p24003, p24006, alcohol, Depression_pmh, p22189, Diabetes_pmh, highest_level, DemFH, PDFH, p26206, p26260, hyperlipidemia_pmh, Hypertension_pmh, obesity, p22032_i0, p31,
                                smoking, p2020_i0,hearing_loss, visual_loss))
missing_summary <- miss_var_summary(Missing)
missing_summary <- as.data.frame(missing_summary)
missing_summary$pct_miss <- as.numeric(missing_summary$pct_miss)
write.csv(missing_summary, "missing_summary_TBI.csv")
columns_to_check <- c("Age_at_recruitment", "alcohol", "Depression_pmh", "Diabetes_pmh", "highest_level", "DemFH", "PDFH", "p26206", "p26260", "hyperlipidemia_pmh", "Hypertension_pmh", "obesity", "p22032_i0", "p31",
                      "smoking", "p2020_i0","hearing_loss", "visual_loss", "p24003", "p24006", "p22189")

M_td_df <- merge(td_df, Adj_1_cov, by = "eid", all.x = T)
F_td_df <- subset(M_td_df, M_td_df$p31=="Female")
M_td_df <- subset(M_td_df, M_td_df$p31=="Male")

crude_incidence <- F_td_df %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )

crude_incidence <- M_td_df %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )



# Create the new variable that indicates missing data
Overall$missing_data_flag <- ifelse(rowSums(is.na(Overall[, columns_to_check])) > 0, 1, 0)
Non<- c("p21022","p22189","p21001_i0", "p24003", "p24006")
myVars <- c("p21022", "p31", "p22189", "highest_level", "smoking", "alcohol","p21000_i0","p22032_i0", "p24003", "p24006",  "p21001_i0","p2178", "Lives_alone", "p2020_i0","DemFH", "PDFH","hearing_loss", "visual_loss" )
catVars <- c("p31", "highest_level","p20116_i0", "p21000_i0","p22032_i0",  "p2178", "Lives_alone", "p2020_i0","hearing_loss", "visual_loss","DemFH", "PDFH")
tab <- CreateTableOne(vars = myVars, data = Overall, strata="missing_data_flag", factorVars = catVars)
tab <- print(tab, formatOptions = list(big.mark = ","), nonnormal =Non, includeNa = T)
tab
write.csv(tab, "Missing_data_demographics_TBI.csv")

Adj_1 <- merge(Adj_1_cov,Simple_3_months_lag, by= "eid", all.y = T)

Adj_2_cov <- Overall %>% select(c(eid, hyperlipidemia_pmh,Hypertension_pmh, Diabetes_pmh, Cerebrovascular_pmh, Cardiovascular_pmh, Depression_pmh, p26206, p26260, DemFH, PDFH, hearing_loss, visual_loss, APOE_ε4_binary))

Adj_2<- merge(Adj_2_cov,Adj_1, by= "eid", all.y = T)

Adj_2 <- Adj_2 %>%
  mutate(Age_at_injury = case_when(
    is.na(Time_to_injury_1) ~ "No Injury",
    Time_to_injury_1 < 40    ~ "Under 40",
    Time_to_injury_1 >= 40 & Time_to_injury_1 < 60 ~ "40 to 60",
    Time_to_injury_1 >= 60   ~ "60+"
  ))

Adj_2_td_df_lag <- tmerge(Adj_2, Adj_2, id = eid, outcome = event(stop, event))
Adj_2_td_df_lag <- tmerge(Adj_2_td_df_lag, Adj_2, id = eid, head_injury = tdc(head_injury_start_1))
Adj_2_td_df_lag$head_injury <- as.factor(Adj_2_td_df_lag$head_injury)

Adj_2_td_df_lag$Age_group <- cut(Adj_2_td_df_lag$Age_at_recruitment, breaks = c(-Inf, 60, Inf), labels = c("<60", ">60"))

convert_variables <- function(df) {
s
  df$head_injury <- factor(df$head_injury)
  df$Age_group <- factor(df$Age_group)
  df$alcohol <- factor(df$alcohol)
  df$Depression_pmh <- factor(df$Depression_pmh)
  df$Diabetes_pmh <- factor(df$Diabetes_pmh)
  df$highest_level <- factor(df$highest_level)
  df$DemFH <- factor(df$DemFH)
  df$PDFH <- factor(df$PDFH)
  df$hyperlipidemia_pmh <- factor(df$hyperlipidemia_pmh)
  df$Hypertension_pmh <- factor(df$Hypertension_pmh)
  df$obesity <- factor(df$obesity)
  df$p22032_i0 <- factor(df$p22032_i0)
  df$p31 <- factor(df$p31)
  df$smoking <- factor(df$smoking)
  df$p2020_i0 <- factor(df$p2020_i0)
  df$visual_loss <- factor(df$visual_loss)
  df$hearing_loss <- factor(df$hearing_loss)
  df$APOE_ε4_binary <- factor(df$APOE_ε4_binary)
  
 
  df$p26206 <- as.numeric(df$p26206)
  df$p26260 <- as.numeric(df$p26260)
  df$p24003 <- as.numeric(df$p24003)
  df$p24006 <- as.numeric(df$p24006)
  
  return(df)
}

Adj_2_td_df_lag <- convert_variables(Adj_2_td_df_lag)
Adj_2_td_df_lag$log_time <- log(Adj_2_td_df_lag$tstop)

com <- subset(Overall, Overall$missing_data_flag==0)
Compl <- subset(Adj_2_td_df_lag, Adj_2_td_df_lag$eid %in% com$eid)
cox_fit_unadj <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = Compl, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_unadj)

cox_fit_B_full <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury + Age_at_recruitment+alcohol+
                          Depression_pmh+p22189+Diabetes_pmh+highest_level+
                          DemFH+PDFH+p26206+p26260+
                          hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Adj_2_td_df_lag, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_B_full)


ph_test <- cox.zph(cox_fit_B_full)
print(ph_test)

export_cox_summary_to_csv(cox_fit_B_full, "cox_fit_non_imputed_full_3_months_TBI.csv")

plot_list <- ggcoxzph(ph_test, plot = FALSE)

for (i in seq_along(plot_list)) {
  plot(plot_list[[i]])
}

arranged_plots <- arrangeGrob(grobs = plot_list, ncol = 1)

ggsave("ggcoxzph.pdf", arranged_plots, width = 8, height = 4 * length(plot_list), limitsize = FALSE)

columns_to_check <- c("Age_at_recruitment","alcohol","Depression_pmh","Diabetes_pmh","highest_level","DemFH","PDFH","p26206","p26260", "p22189",
                      "hyperlipidemia_pmh","Hypertension_pmh","obesity","p22032_i0","p31","smoking","p2020_i0","visual_loss", "p24003", "p24006", "hearing_loss")
complete_rows <- sum(complete.cases(Adj_2[ , columns_to_check]))
print(paste("Number of rows with no missing values in the specified columns:", complete_rows))

filtered_data <- Adj_2_td_df_lag %>%
  filter(    !is.na(APOE_ε4_binary) &
               !is.na(visual_loss) &
               !is.na(hearing_loss) &
               !is.na(p24003) &
               !is.na(p24006) &
               !is.na(p22189) &
               !is.na(Age_at_recruitment) & 
               !is.na(alcohol) & 
               !is.na(Depression_pmh) & 
               !is.na(Diabetes_pmh) & 
               !is.na(highest_level) & 
               !is.na(DemFH) & 
               !is.na(PDFH) & 
               !is.na(p26206) & 
               !is.na(p26260) & 
               !is.na(hyperlipidemia_pmh) & 
               !is.na(Hypertension_pmh) & 
               !is.na(obesity) & 
               !is.na(p22032_i0) & 
               !is.na(p31) & 
               !is.na(smoking) & 
               !is.na(p2020_i0))

table(filtered_data$Age_at_injury)

##With APOE status

cox_fit_B_APOE <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury*APOE_ε4_binary + Age_at_recruitment+alcohol+
                          Depression_pmh+p22189+Diabetes_pmh+highest_level+
                          DemFH+PDFH+p26206+p26260+
                          hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Adj_2_td_df_lag, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_B_APOE)

# Extract estimates from the model summary
coef_summary <- summary(cox_fit_B_APOE)$coefficients
vcov_mat <- vcov(cox_fit_B_APOE)


beta_headinjury <- coef_summary["head_injury1", "coef"]
beta_interaction <- coef_summary["head_injury1:APOE_ε4_binary1", "coef"]


HR_APOE0 <- exp(beta_headinjury)  
HR_APOE1 <- exp(beta_headinjury + beta_interaction)  


cat("Hazard ratio for APOE0 (head injury):", HR_APOE0, "\n")
cat("Hazard ratio for APOE1 (head injury):", HR_APOE1, "\n")


var_APOE1 <- vcov_mat["head_injury1", "head_injury1"] +
  vcov_mat["head_injury1:APOE_ε4_binary1", "head_injury1:APOE_ε4_binary1"] +
  2 * vcov_mat["head_injury1", "head_injury1:APOE_ε4_binary1"]


se_APOE1 <- sqrt(var_APOE1)


se_APOE0 <- coef_summary["head_injury1", "se(coef)"]
CI_APOEO <- exp(beta_headinjury + c(-1.96, 1.96) * se_APOE0)


CI_APOE1 <- exp((beta_headinjury + beta_interaction) + c(-1.96, 1.96) * se_APOE1)


cat("95% CI for APOE0:", round(CI_APOEO[1], 3), "-", round(CI_APOEO[2], 3), "\n")
cat("95% CI for APOE1:", round(se_APOE1[1], 3), "-", round(se_APOE1[2], 3), "\n")


##Age seperated

under_40 <- subset(Adj_2_td_df_lag, Adj_2_td_df_lag$Age_at_injury=="No Injury" | Adj_2_td_df_lag$Age_at_injury=="Under 40")

cox_under_40_full <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury + Age_at_recruitment+alcohol+
                             Depression_pmh+p22189+Diabetes_pmh+highest_level+
                             DemFH+PDFH+p26206+p26260+
                             hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = under_40, id = eid, ties='efron' , timefix = TRUE)
summary(cox_under_40_full)

under_60 <- subset(Adj_2_td_df_lag, Adj_2_td_df_lag$Age_at_injury=="No Injury" | Adj_2_td_df_lag$Age_at_injury=="40 to 60")

cox_under_60_full <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury + Age_at_recruitment+alcohol+
                             Depression_pmh+p22189+Diabetes_pmh+highest_level+
                             DemFH+PDFH+p26206+p26260+
                             hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = under_60, id = eid, ties='efron' , timefix = TRUE)
summary(cox_under_60_full)


over_60 <- subset(Adj_2_td_df_lag, Adj_2_td_df_lag$Age_at_injury=="No Injury" | Adj_2_td_df_lag$Age_at_injury=="60+")

cox_over_60_full <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury + Age_at_recruitment+alcohol+
                            Depression_pmh+p22189+Diabetes_pmh+highest_level+
                            DemFH+PDFH+p26206+p26260+
                            hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = over_60, id = eid, ties='efron' , timefix = TRUE)
summary(cox_over_60_full)


##Age seperated

cox_fit_B_interaction <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury*Age_group+alcohol +
                                 Depression_pmh+p22189+Diabetes_pmh+highest_level+
                                 DemFH+PDFH+p26206+p26260+
                                 hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006,
                               data = Adj_2_td_df_lag, 
                               id = eid, ties='efron', timefix = TRUE)

summary(cox_fit_B_interaction)

coef_summary <- summary(cox_fit_B_interaction)$coefficients
vcov_mat <- vcov(cox_fit_B_interaction)
# Extract estimates
beta_headinjury  <- coef_summary["head_injury1", "coef"]
beta_interaction_60  <- coef_summary["head_injury1:Age_group>60", "coef"]

# HR for reference age group (<55)
HR_under60 <- exp(beta_headinjury)

# HR for age group >55
HR_over60 <- exp(beta_headinjury + beta_interaction_60)

# Print results
cat("HR for <60:", HR_under60, "\n")
cat("HR for >60:", HR_over60, "\n")

# Variance for HR in >55
var_over60 <- vcov_mat["head_injury1", "head_injury1"] + 
  vcov_mat["head_injury1:Age_group>60", "head_injury1:Age_group>60"] + 
  2 * vcov_mat["head_injury1", "head_injury1:Age_group>60"]

# Standard error
se_over60 <- sqrt(var_over60)

# 95% CI for <55
CI_under60 <- exp(beta_headinjury + c(-1.96, 1.96) * coef_summary["head_injury1", "se(coef)"])

# 95% CI for >55
CI_over60 <- exp((beta_headinjury + beta_interaction_60) + c(-1.96, 1.96) * se_over60)

# Print results
cat("95% CI for <60:", CI_under60, "\n")
cat("95% CI for >60:", CI_over60, "\n")


##Sex seperated

cox_fit_B_interaction <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury * p31 + 
                                 Age_at_recruitment + alcohol +
                                 Depression_pmh + p22189 + Diabetes_pmh + highest_level +
                                 DemFH + PDFH + p26206 + p26260 +
                                 hyperlipidemia_pmh + Hypertension_pmh + obesity + p22032_i0  + smoking + p2020_i0 + visual_loss + hearing_loss + 
                                 p24003 + p24006, 
                               data = Adj_2_td_df_lag, 
                               id = eid, ties='efron', timefix = TRUE)


summary(cox_fit_B_interaction)

##
coef_summary <- summary(cox_fit_B_interaction)$coefficients
vcov_mat <- vcov(cox_fit_B_interaction)

##

# Extract estimates from the model summary
coef_summary <- summary(cox_fit_B_interaction)$coefficients
vcov_mat <- vcov(cox_fit_B_interaction)


beta_headinjury <- coef_summary["head_injury1", "coef"]
beta_interaction <- coef_summary["head_injury1:p31Male", "coef"]


HR_female <- exp(beta_headinjury)  # HR for females (reference group)
HR_male <- exp(beta_headinjury + beta_interaction)  # HR for males


cat("Hazard ratio for females (head injury):", HR_female, "\n")
cat("Hazard ratio for males (head injury):", HR_male, "\n")


var_male <- vcov_mat["head_injury1", "head_injury1"] +
  vcov_mat["head_injury1:p31Male", "head_injury1:p31Male"] +
  2 * vcov_mat["head_injury1", "head_injury1:p31Male"]


se_male <- sqrt(var_male)


se_female <- coef_summary["head_injury1", "se(coef)"]
CI_female <- exp(beta_headinjury + c(-1.96, 1.96) * se_female)


CI_male <- exp((beta_headinjury + beta_interaction) + c(-1.96, 1.96) * se_male)


cat("95% CI for females:", round(CI_female[1], 3), "-", round(CI_female[2], 3), "\n")
cat("95% CI for males:", round(CI_male[1], 3), "-", round(CI_male[2], 3), "\n")

data <- data.frame(
  variable = c("No Head Injury", " ", "Females", "Males"),
  HR = c(1, NA, 2.08, 1.62), 
  LL = c(NA, NA, 1.76, 1.41),
  UL = c(NA, NA, 2.47, 1.86)
)

data$variable <- factor(data$variable, levels = c("No Head Injury", " ", "Females", "Males"))

Gender_3_months_TBI_no_multiple <- ggplot(data, aes(x = HR, y = variable)) +
  geom_point(aes(fill = variable), shape = 21, size = 4, color = "black", na.rm = TRUE) +  
  geom_errorbar(aes(xmin = LL, xmax = UL), width = 0.1, na.rm = TRUE) + 
  geom_vline(xintercept = 1, linetype = "dotted", color = "red") +
  coord_flip() +
  labs(x = "Hazard Ratio (95% CI)", y = "", 
       title = "") +
  theme_classic(base_size = 15) +      
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(size = 25, hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 15),  
        legend.position = "none") +
  scale_fill_manual(values = c("Males" = "steelblue", "Females" = "darkorange", "No Head Injury" = "grey60")) +
  scale_y_discrete(labels = c("No Head Injury", "Females", "Males"), breaks = c("No Head Injury", "Females", "Males")) +  
  xlim(0.5, 3)
Gender_3_months_TBI_no_multiple

##Multiple head injuries

##How many

Adj_2 <- Adj_2 %>%
  mutate(head_injury_category_check = case_when(
    is.na(head_injury_start_1) & is.na(head_injury_start_2) ~ "No Head Injury",
    !is.na(head_injury_start_1) & is.na(head_injury_start_2) ~ "One Head Injury",
    !is.na(head_injury_start_1) & !is.na(head_injury_start_2) ~ "2+ Head Injuries"
  ))

table(Adj_2$head_injury_category_check)

Adj_2_td_df_multi <- tmerge(Adj_2, Adj_2, id = eid, outcome = event(stop, event))


Adj_2_td_df_multi <- tmerge(Adj_2_td_df_multi, Adj_2, id = eid, 
                            head_injury_1 = tdc(head_injury_start_1))

Adj_2_td_df_multi <- tmerge(Adj_2_td_df_multi, Adj_2, id = eid, 
                            head_injury_2 = tdc(head_injury_start_2))


Adj_2_td_df_multi <- Adj_2_td_df_multi %>%
  mutate(
    head_injury_status = case_when(
      head_injury_1 == 0 & (is.na(head_injury_2) | head_injury_2 == 0) ~ "No head injury",
      head_injury_1 == 1 & (is.na(head_injury_2) | head_injury_2 == 0) ~ "1 head injury",
      head_injury_1 == 1 & head_injury_2 == 1 ~ "2+ head injuries"
    ),
    head_injury_status = factor(head_injury_status, levels = c("No head injury", "1 head injury", "2+ head injuries"))
  )


Adj_2_td_df_multi$Age_group <- cut(
  Adj_2_td_df_multi$Age_at_recruitment, 
  breaks = c(-Inf, 50, 60, Inf), 
  labels = c("<50", "50-60", ">60")
)



cox_fit_B_Adj_2_effects_multi <- coxph(
  formula = Surv(tstart, tstop, outcome) ~ 
    head_injury_status + Age_at_recruitment+alcohol+
    Depression_pmh+p22189+Diabetes_pmh+highest_level+
    DemFH+PDFH+p26206+p26260+
    hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006,
  data = Adj_2_td_df_multi, 
  id = eid, 
  ties = 'efron', 
  timefix = TRUE
)

# Print model summary
summary(cox_fit_B_Adj_2_effects_multi)

# Function to export Cox model summary to CSV 
export_cox_summary_to_csv(cox_fit_B_Adj_2_effects_multi, "cox_fit_non_imputed_full_multiple_TBI.csv")



##Linear test for trend

Adj_2_td_df_multi$head_injury_status_numeric <- recode(Adj_2_td_df_multi$head_injury_status, 
                                                       "No head injury" = 0, 
                                                       "1 head injury" = 1, 
                                                       "2+ head injuries" = 2)

cox_model_linear_trend <- coxph(
  formula = Surv(tstart, tstop, outcome) ~ 
    head_injury_status_numeric + Age_at_recruitment+alcohol+
    Depression_pmh+p22189+Diabetes_pmh+highest_level+
    DemFH+PDFH+p26206+p26260+
    hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006,
  data = Adj_2_td_df_multi, 
  id = eid, 
  ties = 'efron', 
  timefix = TRUE
)

summary(cox_model_linear_trend)

summary(cox_model_linear_trend)$coef["head_injury_status_numeric", "Pr(>|z|)"]
####Presentation

data <- data.frame(
  variable = c("No Head Injury", "Single Head Injury", "2+ Head Injuries"),
  HR = c(1, 1.71, 2.62),
  LL = c(NA, 1.53, 1.91),
  UL = c(NA, 1.91, 3.59),
  Condition = "Neurodegenerative Disease" 
)

data$variable <- factor(data$variable, levels = c("No Head Injury", "Single Head Injury", "2+ Head Injuries"))

Neur_plot_3_months_TBI <- ggplot(data, aes(x = HR, y = variable, color = Condition)) +
  geom_point(shape = 21, size = 4, fill = "black") +  
  geom_errorbar(aes(xmin = LL, xmax = UL), width = 0.1) + 
  geom_vline(aes(xintercept = 1, color = "Reference"), linetype = "dotted", linewidth = 1) +  # Reference line in legend
  coord_flip() +
  labs(x = "Hazard Ratio (95% CI)", y = "", title = "",
       color = "Neurodegenerative Disease") +  
  scale_color_manual(values = c("Neurodegenerative Disease" = "black")) +  # Match earlier plot colors
  theme_classic(base_size = 15) +  
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(size = 25, hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = c(0.3, 0.8)) + 
  xlim(0.5, 4)

Neur_plot_3_months_TBI

##Time from injury


##sex seperated


cox_fit_B_Adj_2_effects_multi <- coxph(
  formula = Surv(tstart, tstop, outcome) ~ 
    head_injury_status*p31 + Age_at_recruitment+alcohol+
    Depression_pmh+p22189+Diabetes_pmh+highest_level+
    DemFH+PDFH+p26206+p26260+
    hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006,
  data = Adj_2_td_df_multi, 
  id = eid, 
  ties = 'efron', 
  timefix = TRUE
)

linearHypothesis(cox_fit_B_Adj_2_effects_multi, c("head_injury_status1 head injury:p31Male = 0", 
                                                  "head_injury_status2+ head injuries:p31Male = 0"))


summary(cox_fit_B_Adj_2_effects_multi)

lh_test <- linearHypothesis(cox_fit_B_Adj_2_effects_multi, 
                            c("head_injury_status1 head injury:p31Male = 0", 
                              "head_injury_status2+ head injuries:p31Male = 0"))

print(lh_test)


coef_summary <- summary(cox_fit_B_Adj_2_effects_multi)$coefficients
vcov_mat <- vcov(cox_fit_B_Adj_2_effects_multi)


beta_injury1 <- coef_summary["head_injury_status1 head injury", "coef"]   # 1 head injury (female)
beta_injury2 <- coef_summary["head_injury_status2+ head injuries", "coef"] # 2+ head injuries (female)


beta_interaction1 <- coef_summary["head_injury_status1 head injury:p31Male", "coef"]   # 1 head injury male interaction
beta_interaction2 <- coef_summary["head_injury_status2+ head injuries:p31Male", "coef"] # 2+ head injuries male interaction

 (HR)
HR_female_injury1 <- exp(beta_injury1)
HR_male_injury1 <- exp(beta_injury1 + beta_interaction1)

HR_female_injury2 <- exp(beta_injury2)
HR_male_injury2 <- exp(beta_injury2 + beta_interaction2)

# Print HRs
cat("HR for females with 1 head injury:", HR_female_injury1, "\n")
cat("HR for males with 1 head injury:", HR_male_injury1, "\n")
cat("HR for females with 2+ head injuries:", HR_female_injury2, "\n")
cat("HR for males with 2+ head injuries:", HR_male_injury2, "\n")

 with 1 head injury
var_male_injury1 <- vcov_mat["head_injury_status1 head injury", "head_injury_status1 head injury"] +
  vcov_mat["head_injury_status1 head injury:p31Male", "head_injury_status1 head injury:p31Male"] +
  2 * vcov_mat["head_injury_status1 head injury", "head_injury_status1 head injury:p31Male"]

 with 2+ head injuries
var_male_injury2 <- vcov_mat["head_injury_status2+ head injuries", "head_injury_status2+ head injuries"] +
  vcov_mat["head_injury_status2+ head injuries:p31Male", "head_injury_status2+ head injuries:p31Male"] +
  2 * vcov_mat["head_injury_status2+ head injuries", "head_injury_status2+ head injuries:p31Male"]

# Standard errors
se_male_injury1 <- sqrt(var_male_injury1)
se_male_injury2 <- sqrt(var_male_injury2)

 with 1 head injury
se_female_injury1 <- coef_summary["head_injury_status1 head injury", "se(coef)"]
CI_female_injury1 <- exp(beta_injury1 + c(-1.96, 1.96) * se_female_injury1)

 with 1 head injury
CI_male_injury1 <- exp((beta_injury1 + beta_interaction1) + c(-1.96, 1.96) * se_male_injury1)

 with 2+ head injuries
se_female_injury2 <- coef_summary["head_injury_status2+ head injuries", "se(coef)"]
CI_female_injury2 <- exp(beta_injury2 + c(-1.96, 1.96) * se_female_injury2)

 with 2+ head injuries
CI_male_injury2 <- exp((beta_injury2 + beta_interaction2) + c(-1.96, 1.96) * se_male_injury2)


cat("95% CI for females with 1 head injury:", round(CI_female_injury1, 3), "\n")
cat("95% CI for males with 1 head injury:", round(CI_male_injury1, 3), "\n")

cat("95% CI for females with 2+ head injuries:", round(CI_female_injury2, 3), "\n")
cat("95% CI for males with 2+ head injuries:", round(CI_male_injury2, 3), "\n")


data <- data.frame(
  variable = c("No Head Injury", " ", 
               "Females Single Head Injury", "Females 2+ Head Injuries", " ", 
               "Males Single Head Injury", "Males 2+ Head Injuries"),  
  HR = c(1, NA, 
         1.94, 4.08, NA, 
         1.58, 2.05),  # NA for spacer rows
  LL = c(NA, NA, 
         1.62, 2.53, NA, 
         1.37, 1.35),  # Lower 95% CI
  UL = c(NA, NA, 
         2.34, 6.58, NA, 
         1.83, 3.12)  # Upper 95% CI
)

data$variable <- factor(data$variable, 
                        levels = c("No Head Injury", " ", 
                                   "Females Single Head Injury", "Females 2+ Head Injuries", " ", 
                                   "Males Single Head Injury", "Males 2+ Head Injuries"))

print(data)

data$variable <- factor(data$variable, 
                        levels = c("No Head Injury", " ", 
                                   "Females Single Head Injury", "Males Single Head Injury", "", 
                                   "Females 2+ Head Injuries", "Males 2+ Head Injuries"))


# Plot the data
Gender_3_months_TBI_multiple <- ggplot(data, aes(x = HR, y = variable)) +
  geom_point(aes(fill = variable), shape = 21, size = 4, color = "black", na.rm = TRUE) +  
  geom_errorbar(aes(xmin = LL, xmax = UL), width = 0.1, na.rm = TRUE) + 
  geom_vline(xintercept = 1, linetype = "dotted", color = "red") +
  coord_flip() +
  labs(x = "Hazard Ratio (95% CI)", y = "", 
       title = "") +
  theme_classic(base_size = 15) +     
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(size = 25, hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +  
  scale_fill_manual(values = c("Males Single Head Injury" = "steelblue", 
                               "Males 2+ Head Injuries" = "steelblue", 
                               "Females Single Head Injury" = "darkorange", 
                               "Females 2+ Head Injuries" = "darkorange", 
                               "No Head Injury" = "grey60")) +
  scale_y_discrete(labels = c("No Head Injury", 
                              "Females Single Head Injury", "Males Single Head Injury", 
                              "Females 2+ Head Injuries", "Males 2+ Head Injuries"), 
                   breaks = c("No Head Injury", 
                              "Females Single Head Injury", "Males Single Head Injury", 
                              "Females 2+ Head Injuries", "Males 2+ Head Injuries")) +  
  xlim(0.5, 7)  # Adjust x-axis limits

Gender_3_months_TBI_multiple

combined_plot <- plot_grid(Gender_3_months_TBI_no_multiple, Gender_3_months_TBI_multiple, ncol = 1, labels="AUTO", label_size = 22,
                           rel_heights = c(1, 1.5))

ggsave("Combined_Multiple_risk_TBI_gener.tif", combined_plot, width = 5, height = 10)


##1 year lag

td_df_lag <- tmerge(Simple_1_year_lag, Simple_1_year_lag, id = eid, outcome = event(stop, event))
td_df_lag <- tmerge(td_df_lag, Simple_1_year_lag, id = eid, head_injury = tdc(head_injury_start_1))
td_df_lag$head_injury <- as.factor(td_df_lag$head_injury)
cox_fit_A <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = td_df_lag, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_A)

sum(td_df_lag$tstop-td_df_lag$tstart)
sum(Simple_1_year_lag$stop)

export_cox_summary_to_csv(cox_fit_A, "cox_fit_Unadjusted_TBI.csv")

crude_incidence <- td_df_lag %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )

###Adjusted     

Adj_1_lag <- merge(Adj_1_cov,Simple_1_year_lag, by= "eid", all.y = T)

Adj_2_cov <- Overall %>% select(c(eid, hyperlipidemia_pmh,Hypertension_pmh,visual_loss, hearing_loss, Diabetes_pmh, Cerebrovascular_pmh, Cardiovascular_pmh, Depression_pmh, p26206, p26260, DemFH, PDFH))

Adj_2_lag <- merge(Adj_2_cov,Adj_1_lag, by= "eid", all.y = T)

Adj_2_lag_td_df_lag <- tmerge(Adj_2_lag, Adj_2_lag, id = eid, outcome = event(stop, event))
Adj_2_lag_td_df_lag <- tmerge(Adj_2_lag_td_df_lag, Adj_2_lag, id = eid, head_injury = tdc(head_injury_start_1))
Adj_2_lag_td_df_lag$head_injury <- as.factor(Adj_2_lag_td_df_lag$head_injury)

Adj_2_lag_td_df_lag$Age_group <- cut(Adj_2_lag_td_df_lag$Age_at_recruitment, breaks = c(-Inf, 50, 60, Inf), labels = c("<50", "50-60", ">60"))
Adj_2_lag_td_df_lag <- convert_variables(Adj_2_lag_td_df_lag)
cox_fit_B_full <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury + Age_at_recruitment+alcohol+
                          Depression_pmh+p22189+Diabetes_pmh+highest_level+
                          DemFH+PDFH+p26206+p26260+
                          hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Adj_2_lag_td_df_lag, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_B_full)

ph_test <- cox.zph(cox_fit_B_full)
print(ph_test)

export_cox_summary_to_csv(cox_fit_B_full, "cox_fit_non_imputed_full_TBI.csv")

plot_list <- ggcoxzph(ph_test, plot = FALSE)

columns_to_check <- c("p31", "p2020_i0", "p22189", "highest_level", "p24003", "p24006", 
                      "p22032_i0", "obesity", "smoking", "alcohol", "hyperlipidemia_pmh", 
                      "Hypertension_pmh", "Diabetes_pmh", "Cerebrovascular_pmh", 
                      "Cardiovascular_pmh", "Depression_pmh", "p26206", "p26260")
complete_rows <- sum(complete.cases(Adj_2_lag[ , columns_to_check]))
print(paste("Number of rows with no missing values in the specified columns:", complete_rows))

##Multiple head injuries


Adj_2_lag_td_df_lag_multi <- tmerge(Adj_2_lag, Adj_2_lag, id = eid, outcome = event(stop, event))


Adj_2_lag_td_df_lag_multi <- tmerge(Adj_2_lag_td_df_lag_multi, Adj_2_lag, id = eid, 
                                    head_injury_1 = tdc(head_injury_start_1))

Adj_2_lag_td_df_lag_multi <- tmerge(Adj_2_lag_td_df_lag_multi, Adj_2_lag, id = eid, 
                                    head_injury_2 = tdc(head_injury_start_2))


Adj_2_lag_td_df_lag_multi <- Adj_2_lag_td_df_lag_multi %>%
  mutate(
    head_injury_status = case_when(
      head_injury_1 == 0 & (is.na(head_injury_2) | head_injury_2 == 0) ~ "No head injury",
      head_injury_1 == 1 & (is.na(head_injury_2) | head_injury_2 == 0) ~ "1 head injury",
      head_injury_1 == 1 & head_injury_2 == 1 ~ "2+ head injuries"
    ),
    head_injury_status = factor(head_injury_status, levels = c("No head injury", "1 head injury", "2+ head injuries"))
  )


Adj_2_lag_td_df_lag_multi$Age_group <- cut(
  Adj_2_lag_td_df_lag_multi$Age_at_recruitment, 
  breaks = c(-Inf, 50, 60, Inf), 
  labels = c("<50", "50-60", ">60")
)


cox_fit_B_Adj_2_effects_lag_With_IPAQ <- coxph(
  formula = Surv(tstart, tstop, outcome) ~ 
    head_injury_status + Age_at_recruitment+alcohol+
    Depression_pmh+p22189+Diabetes_pmh+highest_level+
    DemFH+PDFH+p26206+p26260+
    hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, 
  data = Adj_2_lag_td_df_lag_multi, 
  id = eid, 
  ties = 'efron', 
  timefix = TRUE
)

# Print model summary
summary(cox_fit_B_Adj_2_effects_lag_With_IPAQ)

# Function to export Cox model summary to CSV 
export_cox_summary_to_csv(cox_fit_B_Adj_2_effects_lag_With_IPAQ, "cox_fit_non_imputed_full_multiple_TBI.csv")


####Presentation

data <- data.frame(
  variable = c("No Head Injury", "Single Head Injury", "2+ Head Injuries"),
  HR = c(1, 1.60, 2.64),
  LL = c(NA, 1.43, 1.97),
  UL = c(NA, 1.78, 3.54)
)
data$variable <- factor(data$variable, levels = c("No Head Injury", "Single Head Injury", "2+ Head Injuries"))
Neur_plot_1_year_TBI <- ggplot(data, aes(x = HR, y = variable)) +
  geom_point(shape = 21, size = 4, fill = "black") +    # Dot for HR
  geom_errorbar(aes(xmin = LL, xmax = UL), width = 0.1) + 
  geom_vline(xintercept = 1, linetype = "dotted", color = "red") +
  coord_flip() +
  labs(x = "Hazard Ratio (95% CI)", y ="", 
       title = "1 year lag") +
  theme_classic(base_size = 15) +       # Minimal theme
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(size = 25, hjust = 0.5),
        axis.title.x = element_text(size = 20))      +xlim(0.5,4)     


Neur_plot_1_year_TBI

###################################3 year lag

td_df_lag_3 <- tmerge(Simple_3_year_lag, Simple_3_year_lag, id = eid, outcome = event(stop, event))
td_df_lag_3 <- tmerge(td_df_lag_3, Simple_3_year_lag, id = eid, head_injury = tdc(head_injury_start_1))
td_df_lag_3$head_injury <- as.factor(td_df_lag_3$head_injury)
cox_fit_B_lag <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = td_df_lag_3, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_B_lag)

export_cox_summary_to_csv(cox_fit_B_lag, "cox_fit_3_year_Unadj_TBI.csv")

crude_incidence <- td_df_lag_3 %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )

#Adjustment set 1: Age_at_recruitment, p31, p2020_i0+visual_loss+hearing_loss+p24003+p24006, p22189, highest_level, p24003, p24006, p22032_i0, obesity, p20117_i0,p20116_i0      

Adj_3_lag <- merge(Adj_1_cov,Simple_3_year_lag, by= "eid", all.y = T)

Adj_3_lag <- merge(Adj_2_cov,Adj_3_lag, by= "eid", all.y = T)

Adj_3_lag_td_df_lag_3 <- tmerge(Adj_3_lag, Adj_3_lag, id = eid, outcome = event(stop, event))
Adj_3_lag_td_df_lag_3 <- tmerge(Adj_3_lag_td_df_lag_3, Adj_3_lag, id = eid, head_injury = tdc(head_injury_start_1))
Adj_3_lag_td_df_lag_3$head_injury <- as.factor(Adj_3_lag_td_df_lag_3$head_injury)

Adj_3_lag_td_df_lag_3$Age_group <- cut(Adj_3_lag_td_df_lag_3$Age_at_recruitment, breaks = c(-Inf, 50, 60, Inf), labels = c("<50", "50-60", ">60"))

Adj_3_lag_td_df_lag_3 <- convert_variables(Adj_3_lag_td_df_lag_3)
cox_fit_B_Adj_3_effects_lag <- coxph(formula = Surv(tstart, tstop, outcome) ~ 
                                       head_injury + Age_at_recruitment+alcohol+
                                       Depression_pmh+p22189+Diabetes_pmh+highest_level+
                                       DemFH+PDFH+p26206+p26260+
                                       hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data =Adj_3_lag_td_df_lag_3, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_B_Adj_3_effects_lag)

ph_test <- cox.zph(cox_fit_B_Adj_3_effects_lag)
print(ph_test)

plot_list <- ggcoxzph(ph_test, plot = FALSE)

# Loop through each plot and print it
for (i in seq_along(plot_list)) {
  plot(plot_list[[i]])
}

arranged_plots <- arrangeGrob(grobs = plot_list, ncol = 1)

ggsave("Three_year_ggcoxzph.pdf", arranged_plots, width = 8, height = 4 * length(plot_list), limitsize = FALSE)

export_cox_summary_to_csv(cox_fit_B_Adj_3_effects_lag, "cox_fit_non_imputed_full_3_year_TBI.csv")

##Multiple head injuries


Adj_3_lag_td_df_lag_multi <- tmerge(Adj_3_lag, Adj_3_lag, id = eid, outcome = event(stop, event))


Adj_3_lag_td_df_lag_multi <- tmerge(Adj_3_lag_td_df_lag_multi, Adj_3_lag, id = eid, 
                                    head_injury_1 = tdc(head_injury_start_1))

Adj_3_lag_td_df_lag_multi <- tmerge(Adj_3_lag_td_df_lag_multi, Adj_3_lag, id = eid, 
                                    head_injury_2 = tdc(head_injury_start_2))

Adj_3_lag_td_df_lag_multi <- Adj_3_lag_td_df_lag_multi %>%
  mutate(
    head_injury_status = case_when(
      head_injury_1 == 0 & (is.na(head_injury_2) | head_injury_2 == 0) ~ "No head injury",
      head_injury_1 == 1 & (is.na(head_injury_2) | head_injury_2 == 0) ~ "1 head injury",
      head_injury_1 == 1 & head_injury_2 == 1 ~ "2+ head injuries"
    ),
    head_injury_status = factor(head_injury_status, levels = c("No head injury", "1 head injury", "2+ head injuries"))
  )

Adj_3_lag_td_df_lag_multi$Age_group <- cut(
  Adj_3_lag_td_df_lag_multi$Age_at_recruitment, 
  breaks = c(-Inf, 50, 60, Inf), 
  labels = c("<50", "50-60", ">60")
)

cox_fit <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury_status, 
                 data = Adj_3_lag_td_df_lag_multi, id = eid, ties = 'efron', timefix = TRUE)

export_cox_summary_to_csv(cox_fit, "cox_fit_3_year_Unadjusted_multiple_TBI.csv")

cox_fit_3_year_multiple <- coxph(
  formula = Surv(tstart, tstop, outcome) ~  head_injury_status + + Age_at_recruitment+alcohol+
    Depression_pmh+p22189+Diabetes_pmh+highest_level+
    DemFH+PDFH+p26206+p26260+
    hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, 
  data = Adj_3_lag_td_df_lag_multi, 
  id = eid, 
  ties = 'efron', 
  timefix = TRUE
)

# Print model summary
summary(cox_fit_3_year_multiple)


# Function to export Cox model summary to CSV 
export_cox_summary_to_csv(cox_fit_3_year_multiple, "cox_fit_non_imputed_full_3_year_multiple_TBI.csv")

data <- data.frame(
  variable = c("No Head Injury", "Single Head Injury", "2+ Head Injuries"),
  HR = c(1, 1.33, 2.05),
  LL = c(NA, 1.18, 1.50),
  UL = c(NA, 1.47, 2.86)
)
data$variable <- factor(data$variable, levels = c("No Head Injury", "Single Head Injury", "2+ Head Injuries"))
Neur_plot_3_year_TBI <- ggplot(data, aes(x = HR, y = variable)) +
  geom_point(shape = 21, size = 4, fill = "black") +    # Dot for HR
  geom_errorbar(aes(xmin = LL, xmax = UL), width = 0.1) + 
  geom_vline(xintercept = 1, linetype = "dotted", color = "red") +
  coord_flip() +
  labs(x = "Hazard Ratio (95% CI)", y ="", 
       title = "3 year lag") +
  theme_classic(base_size = 15) +       # Minimal theme
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(size = 25, hjust = 0.5),
        axis.title.x = element_text(size = 20))      +xlim(0.5,4)     


Neur_plot_3_year_TBI

combined_plot <- plot_grid(Neur_plot_1_year_TBI, Neur_plot_3_year_TBI, ncol = 1)

additional_data <- data.frame(
  variable = c("No Head Injury", "Single Head Injury", "2+ Head Injuries"),
  Parkinsons_HR = c(1,1.17, 2.34),
  Parkinsons_LL = c(NA, 0.93, 1.36),
  Parkinsons_UL = c(NA, 1.47, 4.04),
  Dementia_HR = c(1,1.91, 3.15),
  Dementia_LL = c(NA, 1.70, 2.27),
  Dementia_UL = c(NA, 2.16, 4.37)
)

additional_data$variable <- factor(additional_data$variable, levels = c("No Head Injury", "Single Head Injury", "2+ Head Injuries"))

# Convert data to long format
plot_data <- additional_data %>%
  pivot_longer(cols = starts_with("Parkinsons_") | starts_with("Dementia_"),
               names_to = c("Condition", "Type"),
               names_sep = "_",
               values_to = "Value") %>%
  pivot_wider(names_from = Type, values_from = Value)


plot_data <- plot_data %>%
  mutate(Group = factor(paste(Condition, variable, sep=" - "), 
                        levels = c("Dementia - No Head Injury", "Parkinsons - No Head Injury", "Dementia - Single Head Injury",
                                   "Parkinsons - Single Head Injury", "Dementia - 2+ Head Injuries", "Parkinsons - 2+ Head Injuries")))

Dem_mulit_plot <- ggplot(plot_data, aes(x = Group, y = HR, color = Condition)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +  
  geom_errorbar(aes(ymin = LL, ymax = UL), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "#F8766D") +
  scale_color_manual(values = c("Parkinsons" = "#00BFC4", "Dementia" = "#F8766D")) +
  labs(y = "Hazard Ratio (95% CI)", x = "", 
       title = "3-Year Lag: Parkinson's vs Dementia Risk") +
  theme_classic(base_size = 15) +
  scale_x_discrete(labels = c("", "No Head Injury", "","One Injury", "","2+ Head Injuries")) +
  theme(axis.text.x = element_blank(),
        plot.title = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.position = c(0.2, 0.65)) +
  ylim(0.5, 5) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = 0.1, ymax = Inf, alpha = 0.1, fill = "black") +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = 0.1, ymax = Inf, alpha = 0.1, fill = "black") +
  annotate("rect", xmin = 4.5, xmax = 6.5, ymin = 0.1, ymax = Inf, alpha = 0.1, fill = "black") +
  annotate("text", x = 1.5, y = -1, label = "No Head Injury", size = 5) +
  annotate("text", x = 3.5, y = -1, label = "One Injury", size = 5) +
  annotate("text", x = 5.5, y = -1, label = "2+ Head Injuries", size = 5)
Dem_mulit_plot

ggsave("Combined_Multiple_risk_TBI.tif", combined_plot, width = 7.5, height = 10)
ggsave("Neur_plot_3_month_TBIs.tif", Neur_plot_3_months_TBI, width = 7.5, height = 3.3)

combined_plot <- plot_grid(Dem_mulit_plot, Neur_plot_3_months_TBI, ncol = 1, 
                           rel_heights = c(1, 1.5),
                           labels = "AUTO")
combined_plot
ggsave("Combined_Multi.tif", combined_plot, width = 7.5, height = 10)


##################################No lag


td_df <- tmerge(Simple, Simple, id = eid, outcome = event(stop, event))
td_df <- tmerge(td_df, Simple, id = eid, head_injury = tdc(head_injury_start_1))
td_df$head_injury <- as.factor(td_df$head_injury)
cox_fit_B_lag <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = td_df, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_B_lag)

export_cox_summary_to_csv(cox_fit_B_lag, "cox_fit_3_year_Unadj_TBI.csv")

crude_incidence <- td_df %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )

#Adjustment set 1: Age_at_recruitment, p31, p2020_i0+visual_loss+hearing_loss+p24003+p24006, p22189, highest_level, p24003, p24006, p22032_i0, obesity, p20117_i0,p20116_i0      

Adj <- merge(Adj_1_cov,Simple, by= "eid", all.y = T)

Adj <- merge(Adj_2_cov,Adj, by= "eid", all.y = T)

Adj_lag_td_df <- tmerge(Adj, Adj, id = eid, outcome = event(stop, event))
Adj_lag_td_df <- tmerge(Adj_lag_td_df, Adj, id = eid, head_injury = tdc(head_injury_start_1))
Adj_lag_td_df$head_injury <- as.factor(Adj_lag_td_df$head_injury)

Adj_lag_td_df$Age_group <- cut(Adj_lag_td_df$Age_at_recruitment, breaks = c(-Inf, 50, 60, Inf), labels = c("<50", "50-60", ">60"))

Adj_lag_td_df <- convert_variables(Adj_lag_td_df)
cox_fit_B_Adj <- coxph(formula = Surv(tstart, tstop, outcome) ~ 
                         head_injury + Age_at_recruitment+alcohol+
                         Depression_pmh+p22189+Diabetes_pmh+highest_level+
                         DemFH+PDFH+p26206+p26260+
                         hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data =Adj_lag_td_df, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_B_Adj)


##################################Neurodegeneration subtypes################################################

########Dementia

Outcome_dataset_dementia <- Outcome_dataset
Outcome_dataset_dementia$all_cause_dementia[is.na(Outcome_dataset_dementia$all_cause_dementia)] <- 0
# Calculate the end time considering neurodegeneration and censoring
Outcome_dataset_dementia$end <- ifelse(is.na(Outcome_dataset_dementia$Time_to_all_cause_dementia), 
                                       Outcome_dataset_dementia$Age_at_censoring, 
                                       Outcome_dataset_dementia$Time_to_all_cause_dementia)

# Define the Age at head injury, handling head injuries prior to recruitment
Outcome_dataset_dementia$head_injury_start <- ifelse(is.na(Outcome_dataset_dementia$Time_to_injury_1), NA, Outcome_dataset_dementia$Time_to_injury_1)
Outcome_dataset_dementia$head_injury_start <- ifelse(Outcome_dataset_dementia$head_injury_start < Outcome_dataset_dementia$Age_at_recruitment, 
                                                     (Outcome_dataset_dementia$Age_at_recruitment+0.01), 
                                                     Outcome_dataset_dementia$head_injury_start)

Outcome_dataset_dementia$head_injury_start <- Outcome_dataset_dementia$head_injury_start - Outcome_dataset_dementia$Age_at_recruitment

Outcome_dataset_dementia$start <- 0
Outcome_dataset_dementia$stop <- Outcome_dataset_dementia$end - Outcome_dataset_dementia$Age_at_recruitment
Outcome_dataset_dementia$event <- Outcome_dataset_dementia$all_cause_dementia
Outcome_dataset_dementia$time_in_study <- Outcome_dataset_dementia$stop
Outcome_dataset_dementia$Time_to_injury_1 <- as.numeric(Outcome_dataset_dementia$Time_to_injury_1)
Outcome_dataset_dementia$Time_to_all_cause_dementia <- as.numeric(Outcome_dataset_dementia$Time_to_all_cause_dementia)

Simple_dementia <- Outcome_dataset_dementia %>% select(c(eid, start, stop, event, head_injury_start_1, head_injury_start_2, time_in_study, stop_dementia, event_dementia, time_in_study_dementia, stop_Parkinsons, event_Parkinsons, time_in_study_Parkinsons, GP_only, Parkinsons, MND, all_cause_dementia))

##3 months

Outcome_dataset_dementia_3_month_lag <- Outcome_dataset_dementia
Outcome_dataset_dementia_3_month_lag$Time_to_injury_1 <- as.numeric(Outcome_dataset_dementia_3_month_lag$Time_to_injury_1)
Outcome_dataset_dementia_3_month_lag$Time_to_all_cause_dementia <- as.numeric(Outcome_dataset_dementia_3_month_lag$Time_to_all_cause_dementia)

initial_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_dementia_3_month_lag[[time_col]]))
})

for (i in 1:2) {
  start_col <- paste("head_injury_start_", i, sep = "")
  gg_col <- paste("Time_to_injury_", i, sep = "")
  time_col <- paste("Time_to_injury_", i, sep = "")
  condition <- Outcome_dataset_dementia_3_month_lag[[time_col]] > (Outcome_dataset_dementia_3_month_lag$Time_to_all_cause_dementia - 0.25)
  Outcome_dataset_dementia_3_month_lag[[start_col]][condition] <- NA
  Outcome_dataset_dementia_3_month_lag[[gg_col]][condition] <- NA
}

updated_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_dementia_3_month_lag[[time_col]]))
})

# Compare before and after
na_changes <- updated_na_counts - initial_na_counts
print(na_changes)

Simple_3_month_lag_dementia <- Outcome_dataset_dementia_3_month_lag %>% select(c(eid, start, stop, event, head_injury_start_1, head_injury_start_2, time_in_study, stop_dementia, event_dementia, time_in_study_dementia, stop_Parkinsons, event_Parkinsons, time_in_study_Parkinsons, GP_only, Parkinsons, MND, all_cause_dementia, Time_to_injury_1))

##1 year

Outcome_dataset_dementia_1_year_lag <- Outcome_dataset_dementia
Outcome_dataset_dementia_1_year_lag$Time_to_injury_1 <- as.numeric(Outcome_dataset_dementia_1_year_lag$Time_to_injury_1)
Outcome_dataset_dementia_1_year_lag$Time_to_all_cause_dementia <- as.numeric(Outcome_dataset_dementia_1_year_lag$Time_to_all_cause_dementia)

initial_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_dementia_1_year_lag[[time_col]]))
})

for (i in 1:2) {
  start_col <- paste("head_injury_start_", i, sep = "")
  time_col <- paste("Time_to_injury_", i, sep = "")
  condition <- Outcome_dataset_dementia_1_year_lag[[time_col]] > (Outcome_dataset_dementia_1_year_lag$Time_to_all_cause_dementia - 1)
  Outcome_dataset_dementia_1_year_lag[[start_col]][condition] <- NA
}

updated_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_1_year_lag[[time_col]]))
})

# Compare before and after
na_changes <- updated_na_counts - initial_na_counts
print(na_changes)

Simple_1_year_lag_dementia <- Outcome_dataset_dementia_1_year_lag %>% select(c(eid, start, stop, event, head_injury_start_1, head_injury_start_2, time_in_study, stop_dementia, event_dementia, time_in_study_dementia, stop_Parkinsons, event_Parkinsons, time_in_study_Parkinsons, GP_only, Parkinsons, MND, all_cause_dementia))

###3 months

###Non adjusted

td_df_dementia <- tmerge(Simple_3_month_lag_dementia, Simple_3_month_lag_dementia, id = eid, outcome = event(stop, event))
td_df_dementia <- tmerge(td_df_dementia, Simple_3_month_lag_dementia, id = eid, head_injury = tdc(head_injury_start_1))
td_df_dementia$head_injury <- as.factor(td_df_dementia$head_injury)
cox_fit_A <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = td_df_dementia, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_A)

export_cox_summary_to_csv(cox_fit_A, "cox_fit_dementia_non_imputed_unadjusted_TBI.csv")

crude_incidence <- td_df_dementia %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )



M_td_df <- merge(td_df_dementia, Adj_1_cov, by = "eid", all.x = T)
F_td_df <- subset(M_td_df, M_td_df$p31=="Female")
M_td_df <- subset(M_td_df, M_td_df$p31=="Male")

crude_incidence <- F_td_df %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )

crude_incidence <- M_td_df %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )


###Adjusted for demographics + genetics

Adj_dementia <- merge(Adj_1_cov,Simple_3_month_lag_dementia, by= "eid", all.y = T)
Adj_dementia <- merge(Adj_2_cov,Adj_dementia, by= "eid", all.y = T)

Adj_dementia <- Adj_dementia %>%
  mutate(Age_at_injury = case_when(
    is.na(Time_to_injury_1) ~ "No Injury",
    Time_to_injury_1 < 40    ~ "Under 40",
    Time_to_injury_1 >= 40 & Time_to_injury_1 < 60 ~ "40 to 60",
    Time_to_injury_1 >= 60   ~ "60+"
  ))



Adj_dementia_df <- tmerge(Adj_dementia, Adj_dementia, id = eid, outcome = event(stop, event))
Adj_dementia_df <- tmerge(Adj_dementia_df, Adj_dementia, id = eid, head_injury = tdc(head_injury_start_1))
Adj_dementia_df$head_injury <- as.factor(Adj_dementia_df$head_injury)

Adj_dementia_df$Age_group <- cut(Adj_dementia_df$Age_at_recruitment, breaks = c(-Inf, 60, Inf), labels = c("<60", ">60"))
cox_fit_dementia_effects_lag_Comp <- coxph(formula = Surv(tstart, tstop, outcome) ~  head_injury + Age_at_recruitment+alcohol+
                                             Depression_pmh+p22189+Diabetes_pmh+highest_level+
                                             DemFH+p26206+
                                             hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Adj_dementia_df, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_dementia_effects_lag_Comp)

filtered_data <- Adj_dementia %>%
  filter(    !is.na(APOE_ε4_binary) &
               !is.na(visual_loss) &
               !is.na(hearing_loss) &
               !is.na(p24003) &
               !is.na(p24006) &
               !is.na(p22189) &
               !is.na(Age_at_recruitment) & 
               !is.na(alcohol) & 
               !is.na(Depression_pmh) & 
               !is.na(Diabetes_pmh) & 
               !is.na(highest_level) & 
               !is.na(DemFH) & 
               !is.na(p26206) & 
               !is.na(hyperlipidemia_pmh) & 
               !is.na(Hypertension_pmh) & 
               !is.na(obesity) & 
               !is.na(p22032_i0) & 
               !is.na(p31) & 
               !is.na(smoking) & 
               !is.na(p2020_i0))

table(filtered_data$p31)

export_cox_summary_to_csv(cox_fit_dementia_effects_lag_Comp, "cox_fit_dementia_non_imputed_adjusted_3_months_TBI.csv")

##With APOE status

cox_fit_B_APOE_demen <- coxph(formula = Surv(tstart, tstop, outcome) ~  head_injury + Age_at_recruitment+alcohol+
                                Depression_pmh+p22189+Diabetes_pmh+highest_level+
                                DemFH+p26206+APOE_ε4_binary+
                                hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Adj_dementia_df, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_B_APOE_demen)


###Age groups

under_40 <- subset(Adj_dementia_df, Adj_dementia_df$Age_at_injury=="No Injury" | Adj_dementia_df$Age_at_injury=="Under 40")

cox_under_40_full <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury + Age_at_recruitment+alcohol+
                             Depression_pmh+p22189+Diabetes_pmh+highest_level+
                             DemFH+p26206+
                             hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = under_40, id = eid, ties='efron' , timefix = TRUE)
summary(cox_under_40_full)

under_60 <- subset(Adj_dementia_df, Adj_dementia_df$Age_at_injury=="No Injury" | Adj_dementia_df$Age_at_injury=="40 to 60")

cox_under_60_full <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury + Age_at_recruitment+alcohol+
                             Depression_pmh+p22189+Diabetes_pmh+highest_level+
                             DemFH+p26206+
                             hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = under_60, id = eid, ties='efron' , timefix = TRUE)
summary(cox_under_60_full)


over_60 <- subset(Adj_dementia_df, Adj_dementia_df$Age_at_injury=="No Injury" | Adj_dementia_df$Age_at_injury=="60+")

cox_over_60_full <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury + Age_at_recruitment+alcohol+
                            Depression_pmh+p22189+Diabetes_pmh+highest_level+
                            DemFH+p26206+
                            hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = over_60, id = eid, ties='efron' , timefix = TRUE)
summary(cox_over_60_full)

##Age Interaction

cox_fit_B_interaction <- coxph(formula = Surv(tstart, tstop, outcome) ~  head_injury * Age_group + p31 +alcohol+
                                 Depression_pmh+p22189+Diabetes_pmh+highest_level+
                                 DemFH+p26206+
                                 hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Adj_dementia_df, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_B_interaction)

coef_summary <- summary(cox_fit_B_interaction)$coefficients
vcov_mat <- vcov(cox_fit_B_interaction)

beta_headinjury  <- coef_summary["head_injury1", "coef"]
beta_interaction_60  <- coef_summary["head_injury1:Age_group>60", "coef"]


HR_under60 <- exp(beta_headinjury)

HR_over60 <- exp(beta_headinjury + beta_interaction_60)


cat("HR for <60:", HR_under60, "\n")
cat("HR for >60:", HR_over60, "\n")


var_over60 <- vcov_mat["head_injury1", "head_injury1"] + 
  vcov_mat["head_injury1:Age_group>60", "head_injury1:Age_group>60"] + 
  2 * vcov_mat["head_injury1", "head_injury1:Age_group>60"]

se_over60 <- sqrt(var_over60)

CI_under60 <- exp(beta_headinjury + c(-1.96, 1.96) * coef_summary["head_injury1", "se(coef)"])


CI_over60 <- exp((beta_headinjury + beta_interaction_60) + c(-1.96, 1.96) * se_over60)


cat("95% CI for <60:", CI_under60, "\n")
cat("95% CI for >60:", CI_over60, "\n")





##Sex seperation

cox_fit_B_interaction <- coxph(formula = Surv(tstart, tstop, outcome) ~  head_injury * p31 + Age_at_recruitment+alcohol+
                                 Depression_pmh+p22189+Diabetes_pmh+highest_level+
                                 DemFH+p26206+
                                 hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Adj_dementia_df, id = eid, ties='efron' , timefix = TRUE)


summary(cox_fit_B_interaction)


##
coef_summary <- summary(cox_fit_B_interaction)$coefficients
vcov_mat <- vcov(cox_fit_B_interaction)

##

coef_summary <- summary(cox_fit_B_interaction)$coefficients
vcov_mat <- vcov(cox_fit_B_interaction)


beta_headinjury <- coef_summary["head_injury1", "coef"]
beta_interaction <- coef_summary["head_injury1:p31Male", "coef"]


HR_female <- exp(beta_headinjury)  # HR for females (reference group)
HR_male <- exp(beta_headinjury + beta_interaction)  # HR for males


cat("Hazard ratio for females (head injury):", HR_female, "\n")
cat("Hazard ratio for males (head injury):", HR_male, "\n")


var_male <- vcov_mat["head_injury1", "head_injury1"] +
  vcov_mat["head_injury1:p31Male", "head_injury1:p31Male"] +
  2 * vcov_mat["head_injury1", "head_injury1:p31Male"]


se_male <- sqrt(var_male)


se_female <- coef_summary["head_injury1", "se(coef)"]
CI_female <- exp(beta_headinjury + c(-1.96, 1.96) * se_female)


CI_male <- exp((beta_headinjury + beta_interaction) + c(-1.96, 1.96) * se_male)


cat("95% CI for females:", round(CI_female[1], 3), "-", round(CI_female[2], 3), "\n")
cat("95% CI for males:", round(CI_male[1], 3), "-", round(CI_male[2], 3), "\n")


###Multiple


Adj_2_lag_dementia_multi <- tmerge(Adj_dementia, Adj_dementia, id = eid, outcome = event(stop, event))


Adj_2_lag_dementia_multi <- tmerge(Adj_2_lag_dementia_multi, Adj_dementia, id = eid, 
                                   head_injury_1 = tdc(head_injury_start_1))

Adj_2_lag_dementia_multi <- tmerge(Adj_2_lag_dementia_multi, Adj_dementia, id = eid, 
                                   head_injury_2 = tdc(head_injury_start_2))


Adj_2_lag_dementia_multi <- Adj_2_lag_dementia_multi %>%
  mutate(
    head_injury_status = case_when(
      head_injury_1 == 0 & (is.na(head_injury_2) | head_injury_2 == 0) ~ "No head injury",
      head_injury_1 == 1 & (is.na(head_injury_2) | head_injury_2 == 0) ~ "1 head injury",
      head_injury_1 == 1 & head_injury_2 == 1 ~ "2+ head injuries"
    ),
    head_injury_status = factor(head_injury_status, levels = c("No head injury", "1 head injury", "2+ head injuries"))
  )


Adj_2_lag_dementia_multi$Age_group <- cut(
  Adj_2_lag_dementia_multi$Age_at_recruitment, 
  breaks = c(-Inf, 50, 60, Inf), 
  labels = c("<50", "50-60", ">60")
)


cox_fit_B_Adj_2_effects_dementia <- coxph(
  formula = Surv(tstart, tstop, outcome) ~ 
    head_injury_status + Age_at_recruitment+alcohol+
    Depression_pmh+p22189+Diabetes_pmh+highest_level+
    DemFH+PDFH+p26206+p26260+
    hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, 
  data = Adj_2_lag_dementia_multi, 
  id = eid, 
  ties = 'efron', 
  timefix = TRUE
)

# Print model summary
summary(cox_fit_B_Adj_2_effects_dementia)

##linear

Adj_2_lag_dementia_multi$head_injury_status_numeric <- recode(Adj_2_lag_dementia_multi$head_injury_status, 
                                                              "No head injury" = 0, 
                                                              "1 head injury" = 1, 
                                                              "2+ head injuries" = 2)

cox_fit_B_Adj_2_effects_dementia <- coxph(
  formula = Surv(tstart, tstop, outcome) ~ 
    head_injury_status_numeric + Age_at_recruitment+alcohol+
    Depression_pmh+p22189+Diabetes_pmh+highest_level+
    DemFH+PDFH+p26206+p26260+
    hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, 
  data = Adj_2_lag_dementia_multi, 
  id = eid, 
  ties = 'efron', 
  timefix = TRUE
)

summary(cox_fit_B_Adj_2_effects_dementia)
##1 year

####Cox Proportional Hazards model#############

###Non adjusted

td_df_lag_dementia <- tmerge(Simple_1_year_lag_dementia, Simple_1_year_lag_dementia, id = eid, outcome = event(stop, event))
td_df_lag_dementia <- tmerge(td_df_lag_dementia, Simple_1_year_lag_dementia, id = eid, head_injury = tdc(head_injury_start_1))
td_df_lag_dementia$head_injury <- as.factor(td_df_lag_dementia$head_injury)
cox_fit_A <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = td_df_lag_dementia, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_A)

export_cox_summary_to_csv(cox_fit_A, "cox_fit_dementia_non_imputed_unadjusted_TBI.csv")

crude_incidence <- td_df_lag_dementia %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )

###Adjusted for demographics + genetics

Adj_dementia_lag_2 <- merge(Adj_1_cov,Simple_1_year_lag_dementia, by= "eid", all.y = T)
Adj_dementia_lag_2 <- merge(Adj_2_cov,Adj_dementia_lag_2, by= "eid", all.y = T)

Adj_dementia_df_lag_2 <- tmerge(Adj_dementia_lag_2, Adj_dementia_lag_2, id = eid, outcome = event(stop, event))
Adj_dementia_df_lag_2 <- tmerge(Adj_dementia_df_lag_2, Adj_dementia_lag_2, id = eid, head_injury = tdc(head_injury_start_1))
Adj_dementia_df_lag_2$head_injury <- as.factor(Adj_dementia_df_lag_2$head_injury)

Adj_dementia_df_lag_2$Age_group <- cut(Adj_dementia_df_lag_2$Age_at_recruitment, breaks = c(-Inf, 50, 60, Inf), labels = c("<50", "50-60", ">60"))

Adj_dementia_df_lag_2 <- convert_variables(Adj_dementia_df_lag_2)

cox_fit_dementia_effects_lag_Comp <- coxph(formula = Surv(tstart, tstop, outcome) ~  head_injury + Age_at_recruitment+alcohol+
                                             Depression_pmh+p22189+Diabetes_pmh+highest_level+
                                             DemFH+p26206+
                                             hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0 , data = Adj_dementia_df_lag_2, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_dementia_effects_lag_Comp)

export_cox_summary_to_csv(cox_fit_dementia_effects_lag_Comp, "cox_fit_dementia_non_imputed_adjusted_TBI.csv")


######## 3 year lag ##############
Outcome_dataset_dementia_3_year_lag <- Outcome_dataset_dementia
Outcome_dataset_dementia_3_year_lag$Time_to_injury_1 <- as.numeric(Outcome_dataset_dementia_3_year_lag$Time_to_injury_1)
Outcome_dataset_dementia_3_year_lag$Time_to_all_cause_dementia <- as.numeric(Outcome_dataset_dementia_3_year_lag$Time_to_all_cause_dementia)

initial_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_dementia_3_year_lag[[time_col]]))
})

for (i in 1:2) {
  start_col <- paste("head_injury_start_", i, sep = "")
  time_col <- paste("Time_to_injury_", i, sep = "")
  condition <- Outcome_dataset_dementia_3_year_lag[[time_col]] > (Outcome_dataset_dementia_3_year_lag$Time_to_all_cause_dementia - 3)
  Outcome_dataset_dementia_3_year_lag[[start_col]][condition] <- NA
}

updated_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_dementia_3_year_lag[[time_col]]))
})

# Compare before and after
na_changes <- updated_na_counts - initial_na_counts
print(na_changes)

Simple_3_year_lag_dementia <- Outcome_dataset_dementia_3_year_lag %>% select(c(eid, start, stop, event, head_injury_start_1, head_injury_start_2, time_in_study, stop_dementia, event_dementia, time_in_study_dementia, stop_Parkinsons, event_Parkinsons, time_in_study_Parkinsons, GP_only, Parkinsons, MND, all_cause_dementia))

####Cox Proportional Hazards model#############

###Non adjusted

td_df_lag_dementia_three <- tmerge(Simple_3_year_lag_dementia, Simple_3_year_lag_dementia, id = eid, outcome = event(stop, event))
td_df_lag_dementia_three <- tmerge(td_df_lag_dementia_three, Simple_3_year_lag_dementia, id = eid, head_injury = tdc(head_injury_start_1))
td_df_lag_dementia_three$head_injury <- as.factor(td_df_lag_dementia_three$head_injury)
cox_fit_A <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = td_df_lag_dementia_three, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_A)

export_cox_summary_to_csv(cox_fit_A, "cox_fit_dementia_non_imputed_unadjusted_three_TBI.csv")

crude_incidence <- td_df_lag_dementia_three %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )


###Adjusted for demographics + genetics

Adj_dementia_lag_2_three <- merge(Adj_1_cov,Simple_3_year_lag_dementia, by= "eid", all.y = T)
Adj_dementia_lag_2_three <- merge(Adj_2_cov,Adj_dementia_lag_2_three, by= "eid", all.y = T)

Adj_dementia_df_lag_2_three <- tmerge(Adj_dementia_lag_2_three, Adj_dementia_lag_2_three, id = eid, outcome = event(stop, event))
Adj_dementia_df_lag_2_three <- tmerge(Adj_dementia_df_lag_2_three, Adj_dementia_lag_2_three, id = eid, head_injury = tdc(head_injury_start_1))
Adj_dementia_df_lag_2_three$head_injury <- as.factor(Adj_dementia_df_lag_2_three$head_injury)

Adj_dementia_df_lag_2_three$Age_group <- cut(Adj_dementia_df_lag_2_three$Age_at_recruitment, breaks = c(-Inf, 50, 60, Inf), labels = c("<50", "50-60", ">60"))


Adj_dementia_df_lag_2_three <- convert_variables(Adj_dementia_df_lag_2_three)


cox_fit_dementia_effects_lag_Comp_three <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury + Age_at_recruitment+alcohol+
                                                   Depression_pmh+p22189+Diabetes_pmh+highest_level+
                                                   DemFH+p26206+
                                                   hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006,  data = Adj_dementia_df_lag_2_three, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_dementia_effects_lag_Comp_three)

export_cox_summary_to_csv(cox_fit_dementia_effects_lag_Comp_three, "cox_fit_dementia_non_imputed_adjusted_three_TBI.csv")

##No lag

####Cox Proportional Hazards model#############

###Non adjusted

td_df_lag_dementia <- tmerge(Simple_dementia, Simple_dementia, id = eid, outcome = event(stop, event))
td_df_lag_dementia <- tmerge(td_df_lag_dementia, Simple_dementia, id = eid, head_injury = tdc(head_injury_start_1))
td_df_lag_dementia$head_injury <- as.factor(td_df_lag_dementia$head_injury)
cox_fit_A <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = td_df_lag_dementia, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_A)

export_cox_summary_to_csv(cox_fit_A, "cox_fit_dementia_non_imputed_unadjusted_TBI.csv")

crude_incidence <- td_df_lag_dementia %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )


###Adjusted for demographics + genetics

Adj_dementia <- merge(Adj_1_cov,Simple_dementia, by= "eid", all.y = T)
Adj_dementia <- merge(Adj_2_cov,Adj_dementia, by= "eid", all.y = T)

Adj_dementia_df <- tmerge(Adj_dementia, Adj_dementia, id = eid, outcome = event(stop, event))
Adj_dementia_df <- tmerge(Adj_dementia_df, Adj_dementia, id = eid, head_injury = tdc(head_injury_start_1))
Adj_dementia_df$head_injury <- as.factor(Adj_dementia_df$head_injury)

Adj_dementia_df$Age_group <- cut(Adj_dementia_df$Age_at_recruitment, breaks = c(-Inf, 50, 60, Inf), labels = c("<50", "50-60", ">60"))

Adj_dementia_df <- convert_variables(Adj_dementia_df)

cox_fit_dementia <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury + Age_at_recruitment+alcohol+
                            Depression_pmh+p22189+Diabetes_pmh+highest_level+
                            DemFH+p26206+
                            hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006,  data = Adj_dementia_df, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_dementia)

export_cox_summary_to_csv(cox_fit_dementia, "cox_fit_dementia_TBI.csv")


#############Parkinsons

Outcome_dataset_Parkinson <- Outcome_dataset
Outcome_dataset_Parkinson$Parkinsons[is.na(Outcome_dataset_Parkinson$Parkinsons)] <- 0


Outcome_dataset_Parkinson$end <- ifelse(is.na(Outcome_dataset_Parkinson$Time_to_PD), 
                                        Outcome_dataset_Parkinson$Age_at_censoring, 
                                        Outcome_dataset_Parkinson$Time_to_PD)


Outcome_dataset_Parkinson$head_injury_start <- ifelse(is.na(Outcome_dataset_Parkinson$Time_to_injury_1), NA, Outcome_dataset_Parkinson$Time_to_injury_1)
Outcome_dataset_Parkinson$head_injury_start <- ifelse(Outcome_dataset_Parkinson$head_injury_start < Outcome_dataset_Parkinson$Age_at_recruitment, 
                                                      (Outcome_dataset_Parkinson$Age_at_recruitment+0.01), 
                                                      Outcome_dataset_Parkinson$head_injury_start)

Outcome_dataset_Parkinson$head_injury_start <- Outcome_dataset_Parkinson$head_injury_start - Outcome_dataset_Parkinson$Age_at_recruitment

Outcome_dataset_Parkinson$start <- 0
Outcome_dataset_Parkinson$stop <- Outcome_dataset_Parkinson$end - Outcome_dataset_Parkinson$Age_at_recruitment
Outcome_dataset_Parkinson$event <- Outcome_dataset_Parkinson$Parkinsons 
Outcome_dataset_Parkinson$time_in_study <- Outcome_dataset_Parkinson$stop
Outcome_dataset_Parkinson$Time_to_injury_1 <- as.numeric(Outcome_dataset_Parkinson$Time_to_injury_1)
Head_inj_within_one_PD <- subset(Outcome_dataset_Parkinson, Outcome_dataset_Parkinson$Time_to_PD < (Outcome_dataset_Parkinson$Time_to_injury_1+1))

Simple_Parkinsons <- Outcome_dataset_Parkinson %>% select(c(eid, start, stop, event, head_injury_start_1, head_injury_start_2, time_in_study, stop_dementia, event_dementia, time_in_study_dementia, stop_Parkinsons, event_Parkinsons, time_in_study_Parkinsons, GP_only, Parkinsons, MND, all_cause_dementia))

##3 months 

Outcome_dataset_Parkinson_3_months_lag <- Outcome_dataset_Parkinson
Outcome_dataset_Parkinson_3_months_lag$Time_to_injury_1 <- as.numeric(Outcome_dataset_Parkinson_3_months_lag$Time_to_injury_1)
Outcome_dataset_Parkinson_3_months_lag$Time_to_PD <- as.numeric(Outcome_dataset_Parkinson_3_months_lag$Time_to_PD)

initial_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_Parkinson_3_months_lag[[time_col]]))
})

for (i in 1:2) {
  start_col <- paste("head_injury_start_", i, sep = "")
  gg_col <- paste("Time_to_injury_", i, sep = "")
  time_col <- paste("Time_to_injury_", i, sep = "")
  condition <- Outcome_dataset_Parkinson_3_months_lag[[time_col]] > (Outcome_dataset_Parkinson_3_months_lag$Time_to_PD - 0.25)
  Outcome_dataset_Parkinson_3_months_lag[[start_col]][condition] <- NA
  Outcome_dataset_Parkinson_3_months_lag[[gg_col]][condition] <- NA
}

updated_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_Parkinson_3_months_lag[[time_col]]))
})

# Compare before and after
na_changes <- updated_na_counts - initial_na_counts
print(na_changes)

Simple_3_month_lag_Parkinsons <- Outcome_dataset_Parkinson_3_months_lag %>% select(c(eid, start, stop, event, head_injury_start_1, head_injury_start_2, time_in_study, stop_dementia, event_dementia, time_in_study_dementia, stop_Parkinsons, event_Parkinsons, time_in_study_Parkinsons, GP_only, Parkinsons, MND, all_cause_dementia, Time_to_injury_1))


##1 year
Outcome_dataset_Parkinson_1_year_lag <- Outcome_dataset_Parkinson
Outcome_dataset_Parkinson_1_year_lag$Time_to_injury_1 <- as.numeric(Outcome_dataset_Parkinson_1_year_lag$Time_to_injury_1)
Outcome_dataset_Parkinson_1_year_lag$Time_to_PD <- as.numeric(Outcome_dataset_Parkinson_1_year_lag$Time_to_PD)

initial_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_Parkinson_1_year_lag[[time_col]]))
})

for (i in 1:2) {
  start_col <- paste("head_injury_start_", i, sep = "")
  time_col <- paste("Time_to_injury_", i, sep = "")
  condition <- Outcome_dataset_Parkinson_1_year_lag[[time_col]] > (Outcome_dataset_Parkinson_1_year_lag$Time_to_PD - 1)
  Outcome_dataset_Parkinson_1_year_lag[[start_col]][condition] <- NA
}

updated_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_Parkinson_1_year_lag[[time_col]]))
})

# Compare before and after
na_changes <- updated_na_counts - initial_na_counts
print(na_changes)

Simple_1_year_lag_Parkinsons <- Outcome_dataset_Parkinson_1_year_lag %>% select(c(eid, start, stop, event, head_injury_start_1, head_injury_start_2, time_in_study, stop_dementia, event_dementia, time_in_study_dementia, stop_Parkinsons, event_Parkinsons, time_in_study_Parkinsons, GP_only, Parkinsons, MND, all_cause_dementia))

####################3 months lag

####Cox Proportional Hazards model#############

###Non adjusted

td_df_Parkinsons <- tmerge(Simple_3_month_lag_Parkinsons, Simple_3_month_lag_Parkinsons, id = eid, outcome = event(stop, event))
td_df_Parkinsons <- tmerge(td_df_Parkinsons, Simple_3_month_lag_Parkinsons, id = eid, head_injury = tdc(head_injury_start_1))
td_df_Parkinsons$head_injury <- as.factor(td_df_Parkinsons$head_injury)
cox_fit_A <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = td_df_Parkinsons, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_A)

export_cox_summary_to_csv(cox_fit_A, "cox_fit_parkinsons_non_imputed_unadjusted_TBI.csv")

crude_incidence <- td_df_Parkinsons %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )

M_td_df <- merge(td_df_Parkinsons, Adj_1_cov, by = "eid", all.x = T)
F_td_df <- subset(M_td_df, M_td_df$p31=="Female")
M_td_df <- subset(M_td_df, M_td_df$p31=="Male")

crude_incidence <- F_td_df %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )

crude_incidence <- M_td_df %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )


###Adjusted for demographics + genetics

Adj_Parkisons_2 <- merge(Adj_1_cov,Simple_3_month_lag_Parkinsons, by= "eid", all.y = T)
Adj_Parkisons_2 <- merge(Adj_2_cov,Adj_Parkisons_2, by= "eid", all.y = T)

Adj_Parkisons_2 <- Adj_Parkisons_2 %>%
  mutate(Age_at_injury = case_when(
    is.na(Time_to_injury_1) ~ "No Injury",
    Time_to_injury_1 < 40    ~ "Under 40",
    Time_to_injury_1 >= 40 & Time_to_injury_1 < 60 ~ "40 to 60",
    Time_to_injury_1 >= 60   ~ "60+"
  ))


Adj_parkisons_df <- tmerge(Adj_Parkisons_2, Adj_Parkisons_2, id = eid, outcome = event(stop, event))
Adj_parkisons_df <- tmerge(Adj_parkisons_df, Adj_Parkisons_2, id = eid, head_injury = tdc(head_injury_start_1))
Adj_parkisons_df$head_injury <- as.factor(Adj_parkisons_df$head_injury)


Adj_parkisons_df$Age_group <- cut(Adj_parkisons_df$Age_at_recruitment, breaks = c(-Inf, 60, Inf), labels = c("<60", ">60"))

Adj_parkisons_df <- convert_variables(Adj_parkisons_df)

cox_fit_parkinsons_effects_lag_Comp <- coxph(formula = Surv(tstart, tstop, outcome) ~  head_injury + Age_at_recruitment+alcohol+
                                               Depression_pmh+p22189+Diabetes_pmh +
                                               PDFH+p26260+
                                               p22032_i0+p31+smoking, data = Adj_parkisons_df, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_parkinsons_effects_lag_Comp)

filtered_data <- Adj_parkisons_df %>%
  filter(    !is.na(APOE_ε4_binary) &
               !is.na(Age_at_recruitment) & 
               !is.na(alcohol) & 
               !is.na(Depression_pmh) & 
               !is.na(p22189) &
               !is.na(Diabetes_pmh) & 
               !is.na(PDFH) & 
               !is.na(p26260) & 
               !is.na(p22032_i0) & 
               !is.na(p31) & 
               !is.na(smoking) )

##With APOE

cox_fit_parkinsons_effects_APOE <- coxph(formula = Surv(tstart, tstop, outcome) ~  head_injury + Age_at_recruitment+alcohol+
                                           Depression_pmh+p22189+Diabetes_pmh +
                                           PDFH+p26260+APOE_ε4_binary+
                                           p22032_i0+p31+smoking, data = Adj_parkisons_df, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_parkinsons_effects_APOE)


##Age separated

###Age groups

under_40 <- subset(Adj_parkisons_df, Adj_parkisons_df$Age_at_injury=="No Injury" | Adj_parkisons_df$Age_at_injury=="Under 40")

cox_under_40_full <- coxph(formula = Surv(tstart, tstop, outcome) ~  head_injury + Age_at_recruitment+alcohol+
                             Depression_pmh+p22189+Diabetes_pmh +
                             PDFH+p26260+
                             p22032_i0+p31+smoking, data = under_40, id = eid, ties='efron' , timefix = TRUE)
summary(cox_under_40_full)

under_60 <- subset(Adj_parkisons_df, Adj_parkisons_df$Age_at_injury=="No Injury" | Adj_parkisons_df$Age_at_injury=="40 to 60")

cox_under_60_full <- coxph(formula = Surv(tstart, tstop, outcome) ~  head_injury + Age_at_recruitment+alcohol+
                             Depression_pmh+p22189+Diabetes_pmh +
                             PDFH+p26260+
                             p22032_i0+p31+smoking, data = under_60, id = eid, ties='efron' , timefix = TRUE)
summary(cox_under_60_full)


over_60 <- subset(Adj_parkisons_df, Adj_parkisons_df$Age_at_injury=="No Injury" | Adj_parkisons_df$Age_at_injury=="60+")

cox_over_60_full <- coxph(formula = Surv(tstart, tstop, outcome) ~  head_injury + Age_at_recruitment+alcohol+
                            Depression_pmh+p22189+Diabetes_pmh +
                            PDFH+p26260+
                            p22032_i0+p31+smoking, data = over_60, id = eid, ties='efron' , timefix = TRUE)
summary(cox_over_60_full)

##Age Interaction

cox_fit_B_interaction <- coxph(formula = Surv(tstart, tstop, outcome) ~  head_injury*Age_group + alcohol+
                                 Depression_pmh+p22189+Diabetes_pmh +
                                 PDFH+p26260+
                                 p22032_i0+p31+smoking, data = Adj_parkisons_df, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_B_interaction)

coef_summary <- summary(cox_fit_B_interaction)$coefficients
vcov_mat <- vcov(cox_fit_B_interaction)

beta_headinjury  <- coef_summary["head_injury1", "coef"]
beta_interaction_60  <- coef_summary["head_injury1:Age_group>60", "coef"]


HR_under60 <- exp(beta_headinjury)


HR_over60 <- exp(beta_headinjury + beta_interaction_60)


cat("HR for <60:", HR_under60, "\n")
cat("HR for >60:", HR_over60, "\n")


var_over60 <- vcov_mat["head_injury1", "head_injury1"] + 
  vcov_mat["head_injury1:Age_group>60", "head_injury1:Age_group>60"] + 
  2 * vcov_mat["head_injury1", "head_injury1:Age_group>60"]


se_over60 <- sqrt(var_over60)


CI_under60 <- exp(beta_headinjury + c(-1.96, 1.96) * coef_summary["head_injury1", "se(coef)"])


CI_over60 <- exp((beta_headinjury + beta_interaction_60) + c(-1.96, 1.96) * se_over60)


cat("95% CI for <60:", CI_under60, "\n")
cat("95% CI for >60:", CI_over60, "\n")



##Sex seperated

Mal_Mul <- subset(Adj_parkisons_df, Adj_parkisons_df$p31=="Male")
Fem_Mul <- subset(Adj_parkisons_df, Adj_parkisons_df$p31=="Female")

cox_fit_parkinsons_inte <- coxph(formula = Surv(tstart, tstop, outcome) ~  head_injury * p31 + Age_at_recruitment+alcohol+
                                   Depression_pmh+p22189+Diabetes_pmh +
                                   PDFH+p26260+
                                   p22032_i0+smoking, data = Adj_parkisons_df, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_parkinsons_inte)



summary(cox_fit_parkinsons_inte)


##
coef_summary <- summary(cox_fit_parkinsons_inte)$coefficients
vcov_mat <- vcov(cox_fit_parkinsons_inte)



beta_headinjury <- coef_summary["head_injury1", "coef"]
beta_interaction <- coef_summary["head_injury1:p31Male", "coef"]


HR_female <- exp(beta_headinjury)  # HR for females (reference group)
HR_male <- exp(beta_headinjury + beta_interaction)  # HR for males


cat("Hazard ratio for females (head injury):", HR_female, "\n")
cat("Hazard ratio for males (head injury):", HR_male, "\n")


var_male <- vcov_mat["head_injury1", "head_injury1"] +
  vcov_mat["head_injury1:p31Male", "head_injury1:p31Male"] +
  2 * vcov_mat["head_injury1", "head_injury1:p31Male"]


se_male <- sqrt(var_male)


se_female <- coef_summary["head_injury1", "se(coef)"]
CI_female <- exp(beta_headinjury + c(-1.96, 1.96) * se_female)


CI_male <- exp((beta_headinjury + beta_interaction) + c(-1.96, 1.96) * se_male)


cat("95% CI for females:", round(CI_female[1], 3), "-", round(CI_female[2], 3), "\n")
cat("95% CI for males:", round(CI_male[1], 3), "-", round(CI_male[2], 3), "\n")


##multi



Adj_2_lag_parkinsons_multi <- tmerge(Adj_Parkisons_2, Adj_Parkisons_2, id = eid, outcome = event(stop, event))


Adj_2_lag_parkinsons_multi <- tmerge(Adj_2_lag_parkinsons_multi, Adj_Parkisons_2, id = eid, 
                                     head_injury_1 = tdc(head_injury_start_1))

Adj_2_lag_parkinsons_multi <- tmerge(Adj_2_lag_parkinsons_multi, Adj_Parkisons_2, id = eid, 
                                     head_injury_2 = tdc(head_injury_start_2))


Adj_2_lag_parkinsons_multi <- Adj_2_lag_parkinsons_multi %>%
  mutate(
    head_injury_status = case_when(
      head_injury_1 == 0 & (is.na(head_injury_2) | head_injury_2 == 0) ~ "No head injury",
      head_injury_1 == 1 & (is.na(head_injury_2) | head_injury_2 == 0) ~ "1 head injury",
      head_injury_1 == 1 & head_injury_2 == 1 ~ "2+ head injuries"
    ),
    head_injury_status = factor(head_injury_status, levels = c("No head injury", "1 head injury", "2+ head injuries"))
  )


Adj_2_lag_parkinsons_multi$Age_group <- cut(
  Adj_2_lag_parkinsons_multi$Age_at_recruitment, 
  breaks = c(-Inf, 50, 60, Inf), 
  labels = c("<50", "50-60", ">60")
)


cox_fit_B_Adj_2_effects_parkinsons <- coxph(
  formula = Surv(tstart, tstop, outcome) ~ 
    head_injury_status + Age_at_recruitment+alcohol+
    Depression_pmh+p22189+Diabetes_pmh +
    PDFH+p26260+
    p22032_i0+p31+smoking, 
  data = Adj_2_lag_parkinsons_multi, 
  id = eid, 
  ties = 'efron', 
  timefix = TRUE
)

# Print model summary
summary(cox_fit_B_Adj_2_effects_parkinsons)

##linear

Adj_2_lag_parkinsons_multi$head_injury_status_numeric <- recode(Adj_2_lag_parkinsons_multi$head_injury_status, 
                                                                "No head injury" = 0, 
                                                                "1 head injury" = 1, 
                                                                "2+ head injuries" = 2)

Numeric_PD <- coxph(
  formula = Surv(tstart, tstop, outcome) ~ 
    head_injury_status_numeric + + Age_at_recruitment+alcohol+
    Depression_pmh+p22189+Diabetes_pmh +
    PDFH+p26260+
    p22032_i0+p31+smoking, 
  data = Adj_2_lag_parkinsons_multi, 
  id = eid, 
  ties = 'efron', 
  timefix = TRUE
)

summary(Numeric_PD)

############################################## 1 year lag

###Non adjusted

td_df_lag_Parkinsons <- tmerge(Simple_1_year_lag_Parkinsons, Simple_1_year_lag_Parkinsons, id = eid, outcome = event(stop, event))
td_df_lag_Parkinsons <- tmerge(td_df_lag_Parkinsons, Simple_1_year_lag_Parkinsons, id = eid, head_injury = tdc(head_injury_start_1))
td_df_lag_Parkinsons$head_injury <- as.factor(td_df_lag_Parkinsons$head_injury)
cox_fit_A <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = td_df_lag_Parkinsons, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_A)

export_cox_summary_to_csv(cox_fit_A, "cox_fit_parkinsons_non_imputed_unadjusted_TBI.csv")

crude_incidence <- td_df_lag_Parkinsons %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )

###Adjusted for demographics + genetics

Adj_Parkisons_lag_2 <- merge(Adj_1_cov,Simple_1_year_lag_Parkinsons, by= "eid", all.y = T)
Adj_Parkisons_lag_2 <- merge(Adj_2_cov,Adj_Parkisons_lag_2, by= "eid", all.y = T)

Adj_parkisons_df_lag_2 <- tmerge(Adj_Parkisons_lag_2, Adj_Parkisons_lag_2, id = eid, outcome = event(stop, event))
Adj_parkisons_df_lag_2 <- tmerge(Adj_parkisons_df_lag_2, Adj_Parkisons_lag_2, id = eid, head_injury = tdc(head_injury_start_1))
Adj_parkisons_df_lag_2$head_injury <- as.factor(Adj_parkisons_df_lag_2$head_injury)


Adj_parkisons_df_lag_2$Age_group <- cut(Adj_parkisons_df_lag_2$Age_at_recruitment, breaks = c(-Inf, 50, 60, Inf), labels = c("<50", "50-60", ">60"))

Adj_parkisons_df_lag_2 <- convert_variables(Adj_parkisons_df_lag_2)

cox_fit_parkinsons_effects_lag_Comp <- coxph(formula = Surv(tstart, tstop, outcome) ~  head_injury + Age_at_recruitment+alcohol+
                                               Depression_pmh+p22189+Diabetes_pmh +
                                               PDFH+p26260+
                                               p22032_i0+p31+smoking, data = Adj_parkisons_df_lag_2, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_parkinsons_effects_lag_Comp)

export_cox_summary_to_csv(cox_fit_parkinsons_effects_lag_Comp, "cox_fit_parkinsons_non_imputed_adjusted_TBI.csv")

########### 3 year lag  ##############################################################################################################

Outcome_dataset_Parkinson_3_year_lag <- Outcome_dataset_Parkinson
Outcome_dataset_Parkinson_3_year_lag$Time_to_injury_1 <- as.numeric(Outcome_dataset_Parkinson_3_year_lag$Time_to_injury_1)
Outcome_dataset_Parkinson_3_year_lag$Time_to_PD <- as.numeric(Outcome_dataset_Parkinson_3_year_lag$Time_to_PD)

initial_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_Parkinson_3_year_lag[[time_col]]))
})

for (i in 1:2) {
  start_col <- paste("head_injury_start_", i, sep = "")
  time_col <- paste("Time_to_injury_", i, sep = "")
  condition <- Outcome_dataset_Parkinson_3_year_lag[[time_col]] > (Outcome_dataset_Parkinson_3_year_lag$Time_to_PD - 3)
  Outcome_dataset_Parkinson_3_year_lag[[start_col]][condition] <- NA
}

updated_na_counts <- sapply(1:2, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Outcome_dataset_Parkinson_3_year_lag[[time_col]]))
})

# Compare before and after
na_changes <- updated_na_counts - initial_na_counts
print(na_changes)

Simple_3_year_lag_Parkinsons <- Outcome_dataset_Parkinson_3_year_lag %>% select(c(eid, start, stop, event, head_injury_start_1, head_injury_start_2, time_in_study, stop_dementia, event_dementia, time_in_study_dementia, stop_Parkinsons, event_Parkinsons, time_in_study_Parkinsons, GP_only, Parkinsons, MND, all_cause_dementia))

####Cox Proportional Hazards model#############

###Non adjusted

td_df_lag_Parkinsons_three <- tmerge(Simple_3_year_lag_Parkinsons, Simple_3_year_lag_Parkinsons, id = eid, outcome = event(stop, event))
td_df_lag_Parkinsons_three <- tmerge(td_df_lag_Parkinsons_three, Simple_3_year_lag_Parkinsons, id = eid, head_injury = tdc(head_injury_start_1))
td_df_lag_Parkinsons_three$head_injury <- as.factor(td_df_lag_Parkinsons_three$head_injury)
cox_fit_A <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = td_df_lag_Parkinsons_three, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_A)

export_cox_summary_to_csv(cox_fit_A, "cox_fit_parkinsons_non_imputed_unadjusted_three_TBI.csv")

crude_incidence <- td_df_lag_Parkinsons_three %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )

###Adjusted for demographics + genetics

Adj_Parkisons_lag_3 <- merge(Adj_1_cov,Simple_3_year_lag_Parkinsons, by= "eid", all.y = T)
Adj_Parkisons_lag_3 <- merge(Adj_2_cov,Adj_Parkisons_lag_3, by= "eid", all.y = T)

Adj_parkisons_df_lag_3 <- tmerge(Adj_Parkisons_lag_3, Adj_Parkisons_lag_3, id = eid, outcome = event(stop, event))
Adj_parkisons_df_lag_3 <- tmerge(Adj_parkisons_df_lag_3, Adj_Parkisons_lag_3, id = eid, head_injury = tdc(head_injury_start_1))
Adj_parkisons_df_lag_3$head_injury <- as.factor(Adj_parkisons_df_lag_3$head_injury)


Adj_parkisons_df_lag_3$Age_group <- cut(Adj_parkisons_df_lag_3$Age_at_recruitment, breaks = c(-Inf, 50, 60, Inf), labels = c("<50", "50-60", ">60"))

Adj_parkisons_df_lag_3 <- convert_variables(Adj_parkisons_df_lag_3)
cox_fit_parkinsons_effects_lag_Comp_three <- coxph(formula = Surv(tstart, tstop, outcome) ~  head_injury + Age_at_recruitment+alcohol+
                                                     Depression_pmh+p22189+Diabetes_pmh +
                                                     PDFH+p26260+
                                                     p22032_i0+p31+smoking, data = Adj_parkisons_df_lag_3, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_parkinsons_effects_lag_Comp_three)

export_cox_summary_to_csv(cox_fit_parkinsons_effects_lag_Comp_three, "cox_fit_parkinsons_non_imputed_adjusted_three_TBI.csv")

#######No lag

####Cox Proportional Hazards model#############

###Non adjusted

td_df_lag_Parkinsons <- tmerge(Simple_Parkinsons, Simple_Parkinsons, id = eid, outcome = event(stop, event))
td_df_lag_Parkinsons <- tmerge(td_df_lag_Parkinsons, Simple_Parkinsons, id = eid, head_injury = tdc(head_injury_start_1))
td_df_lag_Parkinsons$head_injury <- as.factor(td_df_lag_Parkinsons$head_injury)
cox_fit_A <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = td_df_lag_Parkinsons, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_A)

export_cox_summary_to_csv(cox_fit_A, "cox_fit_parkinsons_non_imputed_unadjusted_TBI.csv")

crude_incidence <- td_df_lag_Parkinsons %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )

###Adjusted for demographics + genetics

Adj_Parkisons_lag <- merge(Adj_1_cov,Simple_Parkinsons, by= "eid", all.y = T)
Adj_Parkisons_lag <- merge(Adj_2_cov,Adj_Parkisons_lag, by= "eid", all.y = T)

Adj_parkisons_df <- tmerge(Adj_Parkisons_lag, Adj_Parkisons_lag, id = eid, outcome = event(stop, event))
Adj_parkisons_df <- tmerge(Adj_parkisons_df, Adj_Parkisons_lag, id = eid, head_injury = tdc(head_injury_start_1))
Adj_parkisons_df$head_injury <- as.factor(Adj_parkisons_df$head_injury)


Adj_parkisons_df$Age_group <- cut(Adj_parkisons_df$Age_at_recruitment, breaks = c(-Inf, 50, 60, Inf), labels = c("<50", "50-60", ">60"))
Adj_parkisons_df <- convert_variables(Adj_parkisons_df)

cox_fit_parkinsons <- coxph(formula = Surv(tstart, tstop, outcome) ~  head_injury + Age_at_recruitment+alcohol+
                              Depression_pmh+p22189+Diabetes_pmh +
                              PDFH+p26260+
                              p22032_i0+p31+smoking, data = Adj_parkisons_df, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_parkinsons)

export_cox_summary_to_csv(cox_fit_parkinsons, "cox_fit_parkinsons_non_imputed_adjusted.csv")

##################Plot the age seperated graph

data <- data.frame(
  AgeCategory = c("No Head Injury", "Under 40", "Under 40", "Under 40", 
                  "No Head Injury", "40 to 60", "40 to 60", "40 to 60", 
                  "No Head Injury", "60 and over", "60 and over", "60 and over"),
  Disease = c("Neurodegenerative disease", "Neurodegenerative disease", "Dementia", "Parkinson's disease",
              "Dementia", "Neurodegenerative disease", "Dementia", "Parkinson's disease",
              "Parkinson's disease", "Neurodegenerative disease", "Dementia", "Parkinson's disease"),
  AdjustedHR = c(1, 1.33, 1.38, 1.16, 
                 1, 1.64, 1.89, 1.11, 
                 1, 2.03, 2.24, 1.42),
  LL = c(NA, 1.02, 1.02, 0.74, 
         NA, 1.34, 1.52, 0.74, 
         NA, 1.77, 1.86, 1.05),
  UL = c(NA, 1.73, 1.86, 1.82, 
         NA, 2.00, 2.35, 1.65, 
         NA, 2.34, 2.69, 1.92)
)

data$AgeCategory <- factor(data$AgeCategory, levels = c("No Head Injury", "Under 40", "40 to 60", "60 and over"))
data$Disease <- factor(data$Disease, levels = c("Neurodegenerative disease", "Parkinson's disease", "Dementia"))

age_plot <- ggplot(data, aes(x = Disease, y = AdjustedHR, color = AgeCategory)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = LL, ymax = UL), width = 0.2, linewidth = 0.8, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, color = 'red', linetype = 'dashed') +
  scale_color_manual(values = c("No Head Injury" = "pink", 
                                "Under 40" = "#4F81BD",  
                                "40 to 60" = "#F5A623", 
                                "60 and over" = "#8B9B3D")) +  
  labs(title = "", 
       x = "Disease", 
       y = "Adjusted Hazard Ratio") +
  ylim(0, 4) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        plot.title = element_blank(),
        axis.title.y = element_text(size = 16),
        #axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = c("right"))+
  scale_x_discrete(labels = c("Neurodegenerative disease" = "Neurodegenerative\n disease", 
                              "Parkinson's disease" = "Parkinson's disease", 
                              "Dementia" = "Dementia"))

ggsave("Age_HR.tif", age_plot, width = 7.5, height = 6)


###Is the time between injury and neurodegeneration different?

imp <- Overall %>% select(c("eid","Time_to_Neurodegeneration", "Time_to_all_cause_dementia","Time_to_PD","Time_to_MND","neurodegeneration","p21022", "Age_at_recruitment"))

imp <- merge(imp,Simple_3_months_lag, by="eid")

imp <- imp %>%
  mutate(Age_at_injury = case_when(
    is.na(Time_to_injury_1) ~ "No Injury",
    Time_to_injury_1 < 40    ~ "Under 40",
    Time_to_injury_1 >= 40 & Time_to_injury_1 < 60 ~ "40 to 60",
    Time_to_injury_1 >= 60   ~ "60+"
  ))

imp$time_inj_to_neuro <- imp$Time_to_Neurodegeneration - imp$Time_to_injury_1
imp$time_inj_to_dem <- imp$Time_to_all_cause_dementia - imp$Time_to_injury_1
imp$time_inj_to_PD <- imp$Time_to_PD - imp$Time_to_injury_1

Non<- c("p21022", "Time_to_all_cause_dementia", "Time_to_Neurodegeneration", "Time_to_PD","Time_to_MND", "time_inj_to_neuro", "Time_to_all_cause_dementia", "Time_to_PD", "time_inj_to_dem", "time_inj_to_PD", "Time_to_injury_1")
myVars <- c("p21022","neurodegeneration", "Time_to_Neurodegeneration","all_cause_dementia","Time_to_all_cause_dementia", "Parkinsons", "Time_to_PD", "MND", "Time_to_MND", "time_inj_to_neuro", "time_inj_to_dem", "time_inj_to_PD", "Time_to_injury_1")
catVars <- c("neurodegeneration","all_cause_dementia", "Parkinsons","MND")

tab <- CreateTableOne(vars = myVars, data = imp, strata= "Age_at_injury", factorVars = catVars)
tab <- print(tab, formatOptions = list(big.mark = ","), nonnormal =Non, includeNa = T)
tab

write.csv(tab, "Neurodegeneration_type_TBI.csv")

write.csv(tab, "Time_to_degen_Age.csv")

#######################################################################Sensitivity and Severity############
########Death as competing risk########################################

##remove those who have died

dead_dataset <- subset(Death_reg, !is.na(Death_reg$p40007_i0))

Dead_lag <- Adj_2[!Adj_2$eid %in% dead_dataset$eid, ]

td_df_lag_DEAD <- tmerge(Dead_lag, Dead_lag, id = eid, outcome = event(stop, event))
td_df_lag_DEAD <- tmerge(td_df_lag_DEAD, Dead_lag, id = eid, head_injury = tdc(head_injury_start_1))
td_df_lag_DEAD$head_injury <- as.factor(td_df_lag_DEAD$head_injury)
cox_fit_B_DEAD <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = td_df_lag_DEAD, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_B_DEAD)

export_cox_summary_to_csv(cox_fit_B_DEAD, "cox_fit_GP_removed_Unadj_TBI.csv")

crude_incidence <- td_df_lag_DEAD %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )


Dead_lag_lag <- tmerge(Dead_lag, Dead_lag, id = eid, outcome = event(stop, event))
Dead_lag_lag <- tmerge(Dead_lag_lag, Dead_lag, id = eid, head_injury = tdc(head_injury_start_1))
Dead_lag_lag$head_injury <- as.factor(Dead_lag_lag$head_injury)

Dead_lag_lag$Age_group <- cut(Dead_lag_lag$Age_at_recruitment, breaks = c(-Inf, 50, 60, Inf), labels = c("<50", "50-60", ">60"))
Dead_lag_lag <- convert_variables(Dead_lag_lag)
cox_fit_B_full <- coxph(formula = Surv(tstart, tstop, outcome) ~  head_injury + Age_at_recruitment+alcohol+
                          Depression_pmh+p22189+Diabetes_pmh+highest_level+
                          DemFH+PDFH+p26206+p26260+
                          hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Dead_lag_lag, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_B_full)

export_cox_summary_to_csv(cox_fit_B_full, "cox_fit_non_imputed_death_sensitivity_TBI.csv")

filtered_data <- Dead_lag %>%
  filter(  !is.na(visual_loss) &
             !is.na(hearing_loss) &
             !is.na(p24003) &
             !is.na(p24006) &
             !is.na(p22189) &
             !is.na(Age_at_recruitment) & 
             !is.na(alcohol) & 
             !is.na(Depression_pmh) & 
             !is.na(Diabetes_pmh) & 
             !is.na(highest_level) & 
             !is.na(DemFH) & 
             !is.na(PDFH) & 
             !is.na(p26206) & 
             !is.na(p26260) & 
             !is.na(hyperlipidemia_pmh) & 
             !is.na(Hypertension_pmh) & 
             !is.na(obesity) & 
             !is.na(p22032_i0) & 
             !is.na(p31) & 
             !is.na(smoking) & 
             !is.na(p2020_i0))



##############################################################################################HES only

library(tidyverse)
library(survival)


Head_injured_HES <- read.csv("~/Daniel_Whitehouse/Head_injured_cohort/output/complete_head_injured_Core_HES.csv")

key <- Adj_2 %>% select(c(eid, Age_at_recruitment, p31, p2020_i0, visual_loss, 
                          hearing_loss, p22189, highest_level, p24003, p24006, 
                          p22032_i0, obesity, smoking, alcohol, DemFH, PDFH, 
                          hyperlipidemia_pmh, Hypertension_pmh, Diabetes_pmh, 
                          Cerebrovascular_pmh, Cardiovascular_pmh, Depression_pmh, 
                          p26206, p26260, start, stop, event, stop_dementia, 
                          event_dementia, time_in_study_dementia, stop_Parkinsons, 
                          event_Parkinsons, time_in_study_Parkinsons))

Key_neuro <- Outcome_dataset %>% select(eid, Time_to_Neurodegeneration, p40007_i0)

Head_injured_HES <- merge(Head_injured_HES, key,       by = "eid", all.y = TRUE)
Head_injured_HES <- merge(Head_injured_HES, Key_neuro, by = "eid", all.y = TRUE)

for (i in 1:7) {
  Head_injured_HES[[paste0("Time_to_injury_", i)]] <- as.numeric(Head_injured_HES[[paste0("Time_to_injury_", i)]])
}
Head_injured_HES$Time_to_Neurodegeneration <- as.numeric(Head_injured_HES$Time_to_Neurodegeneration)
Head_injured_HES$p40007_i0                 <- as.numeric(Head_injured_HES$p40007_i0)


Head_injured_HES$head_injured_death <- 0L

for (i in 1:7) {
  time_col <- paste0("Time_to_injury_", i)
  Head_injured_HES$head_injured_death <- ifelse(
    !is.na(Head_injured_HES[[time_col]]) &
      !is.na(Head_injured_HES$p40007_i0) &
      abs(Head_injured_HES[[time_col]] - Head_injured_HES$p40007_i0) <= 0.25,
    1L,
    Head_injured_HES$head_injured_death
  )
}


initial_na_counts <- numeric(7)
updated_na_counts <- numeric(7)

for (i in 1:7) {
  time_col  <- paste0("Time_to_injury_", i)
  start_col <- paste0("head_injury_start_", i)
  sev_col   <- paste0("severity_", i)
  

  Head_injured_HES[[start_col]] <- ifelse(
    is.na(Head_injured_HES[[time_col]]),
    NA,
    Head_injured_HES[[time_col]]
  )
  
  
  Head_injured_HES[[start_col]] <- ifelse(
    !is.na(Head_injured_HES[[time_col]]) &
      !is.na(Head_injured_HES$p40007_i0) &
      abs(Head_injured_HES[[time_col]] - Head_injured_HES$p40007_i0) <= 0.25,
    NA,
    Head_injured_HES[[start_col]]
  )
  
  
  Head_injured_HES[[start_col]] <- ifelse(
    !is.na(Head_injured_HES[[start_col]]) &
      Head_injured_HES[[start_col]] < Head_injured_HES$Age_at_recruitment,
    Head_injured_HES$Age_at_recruitment + 0.01,
    Head_injured_HES[[start_col]]
  )
  
 
  Head_injured_HES[[start_col]] <- Head_injured_HES[[start_col]] - Head_injured_HES$Age_at_recruitment
  

  initial_na_counts[i] <- sum(is.na(Head_injured_HES[[start_col]]))
  
  
  condition <- !is.na(Head_injured_HES[[time_col]]) &
    !is.na(Head_injured_HES$Time_to_Neurodegeneration) &
    Head_injured_HES[[time_col]] > (Head_injured_HES$Time_to_Neurodegeneration - 0.25)
  
  Head_injured_HES[[start_col]][condition] <- NA
  Head_injured_HES[[sev_col]][condition]   <- NA
  

  updated_na_counts[i] <- sum(is.na(Head_injured_HES[[start_col]]))
}

na_changes <- updated_na_counts - initial_na_counts
names(na_changes) <- paste0("head_injury_start_", 1:7)
print(na_changes)

Simple_3_month_HES <- Head_injured_HES %>% 
  select(eid, start, stop, event,
         stop_dementia, event_dementia, time_in_study_dementia,
         stop_Parkinsons, event_Parkinsons, time_in_study_Parkinsons,
         paste0("head_injury_start_", 1:7),
         paste0("severity_", 1:7),
         head_injured_death,
         Age_at_recruitment, p31, p2020_i0, visual_loss, hearing_loss,
         p22189, highest_level, p24003, p24006, p22032_i0,
         obesity, smoking, alcohol, DemFH, PDFH,
         hyperlipidemia_pmh, Hypertension_pmh, Diabetes_pmh,
         Cerebrovascular_pmh, Cardiovascular_pmh, Depression_pmh,
         p26206, p26260)

for (i in 1:7) {
  start_col <- paste0("head_injury_start_", i)
  sev_col   <- paste0("severity_", i)
  Simple_3_month_HES[[sev_col]] <- if_else(
    is.na(Simple_3_month_HES[[start_col]]),
    NA_character_,
    Simple_3_month_HES[[sev_col]]
  )
}


HES_check     <- subset(Simple_3_month_HES,  !is.na(head_injury_start_1))
Overall_check <- subset(Simple_3_months_lag, !is.na(head_injury_start_1))
Comm_on_check <- subset(Overall_check, !Overall_check$eid %in% HES_check$eid)

Simple_HES_sev <- Simple_3_month_HES %>%
  mutate(
    head_injury_start_1 = if_else(
      is.na(head_injury_start_1) & eid %in% Comm_on_check$eid,
      Comm_on_check$head_injury_start_1[match(eid, Comm_on_check$eid)],
      head_injury_start_1
    ),
    severity_1 = if_else(
      is.na(severity_1) & eid %in% Comm_on_check$eid,
      "community",
      severity_1
    )
  )

cat("Community cases added:     ", nrow(Comm_on_check), "\n")
cat("Non-NA injury 1 (before):  ", sum(!is.na(Simple_3_month_HES$head_injury_start_1)), "\n")
cat("Non-NA injury 1 (after):   ", sum(!is.na(Simple_HES_sev$head_injury_start_1)), "\n")


injuries_long <- Simple_HES_sev %>%
  select(eid,
         paste0("head_injury_start_", 1:7),
         paste0("severity_", 1:7)) %>%
  pivot_longer(
    cols          = -eid,
    names_to      = c(".value", "injury_num"),
    names_pattern = "^(.+)_(\\d+)$"
  ) %>%
  rename(injury_time = head_injury_start) %>%
  filter(!is.na(injury_time)) %>%
  arrange(eid, injury_time)


table(injuries_long$severity, useNA = "always")

base <- Simple_HES_sev %>%
  select(eid, start, stop, event,
         stop_dementia, event_dementia,
         stop_Parkinsons, event_Parkinsons,
         head_injured_death,
         Age_at_recruitment, p31, p2020_i0, visual_loss, hearing_loss,
         p22189, highest_level, p24003, p24006, p22032_i0,
         obesity, smoking, alcohol, DemFH, PDFH,
         hyperlipidemia_pmh, Hypertension_pmh, Diabetes_pmh,
         Cerebrovascular_pmh, Cardiovascular_pmh, Depression_pmh,
         p26206, p26260)

td_df_HES <- tmerge(base, base, id = eid,
                    outcome = event(stop, event))

td_df_HES <- tmerge(td_df_HES, injuries_long, id = eid,
                    severity_tv = tdc(injury_time, severity))

# Fill NA (no injury / pre-injury) as "no_injury"
td_df_HES$severity_tv <- if_else(is.na(td_df_HES$severity_tv), "no_injury", td_df_HES$severity_tv)

# Adjust levels to match your table() output above
td_df_HES$severity_tv <- factor(td_df_HES$severity_tv,
                                levels = c("no_injury", "community", "mild", "moderate/severe"))

table(td_df_HES$severity_tv, useNA = "always")  # confirm no NAs

cox_unadj_HES <- coxph(
  Surv(tstart, tstop, outcome) ~ severity_tv,
  data = td_df_HES, id = eid, ties = "efron", timefix = TRUE
)
summary(cox_unadj_HES)

cox_adj_HES <- coxph(
  Surv(tstart, tstop, outcome) ~ severity_tv 
    + Age_at_recruitment+alcohol+
    Depression_pmh+p22189+Diabetes_pmh+highest_level+
    DemFH+PDFH+p26206+p26260+
    hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006,
  data = td_df_HES, id = eid, ties = "efron", timefix = TRUE
)
summary(cox_adj_HES)


Simple_HES_sev %>%
  mutate(severity_group = case_when(
    !is.na(severity_1) ~ severity_1,
    TRUE               ~ "no_injury"
  )) %>%
  group_by(severity_group) %>%
  summarise(
    n          = n(),
    n_outcome  = sum(event == 1, na.rm = TRUE),
    pct_outcome = round(n_outcome / n * 100, 1)
  ) %>%
  arrange(factor(severity_group, levels = c("no_injury", "community", "mild", "moderate/severe")))

td_df_HES %>%
  group_by(eid) %>%
  filter(n() >= 3) %>%
  slice_head(n = 3) %>%  # show first 3 rows per person
  select(eid, tstart, tstop, severity_tv, outcome) %>%
  head(30) 


#################################Only those with GP#########################################################

GP_pres <- read.csv("~/Daniel_Whitehouse/Neurodegeneration/input/GP_record_available.csv")
GP_pres <- subset(GP_pres, !is.na(GP_pres$p42040))
GP_pos <- subset(Adj_2, Adj_2$eid %in% GP_pres$eid)

td_df_lag_GP_pres <- tmerge(GP_pos, GP_pos, id = eid, outcome = event(stop, event))
td_df_lag_GP_pres <- tmerge(td_df_lag_GP_pres, GP_pos, id = eid, head_injury = tdc(head_injury_start_1))
td_df_lag_GP_pres$head_injury <- as.factor(td_df_lag_GP_pres$head_injury)
cox_fit_A <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = td_df_lag_GP_pres, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_A)

export_cox_summary_to_csv(cox_fit_A, "cox_fit_Unadjusted_with_GP_TBI.csv")

crude_incidence <- td_df_lag_GP_pres %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )

Adj_GPPres_lag_td_df_2 <- tmerge(GP_pos, GP_pos, id = eid, outcome = event(stop, event))
Adj_GPPres_lag_td_df_2 <- tmerge(Adj_GPPres_lag_td_df_2, GP_pos, id = eid, head_injury = tdc(head_injury_start_1))
Adj_GPPres_lag_td_df_2$head_injury <- as.factor(Adj_GPPres_lag_td_df_2$head_injury)


Adj_GPPres_lag_td_df_2$Age_group <- cut(Adj_GPPres_lag_td_df_2$Age_at_recruitment, breaks = c(-Inf, 50, 60, Inf), labels = c("<50", "50-60", ">60"))
cox_fit_GPpres_effects_lag_Comp <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury + Age_at_recruitment+alcohol+
                                           Depression_pmh+p22189+Diabetes_pmh+highest_level+
                                           DemFH+PDFH+p26206+p26260+
                                           hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Adj_GPPres_lag_td_df_2, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_GPpres_effects_lag_Comp)

export_cox_summary_to_csv(cox_fit_GPpres_effects_lag_Comp, "cox_fit_Ajd_GPpresent_adjusted_TBI.csv")

filtered_data <- GP_pos %>%
  filter(  !is.na(visual_loss) &
             !is.na(hearing_loss) &
             !is.na(p24003) &
             !is.na(p24006) &
             !is.na(p22189) &
             !is.na(Age_at_recruitment) & 
             !is.na(alcohol) & 
             !is.na(Depression_pmh) & 
             !is.na(Diabetes_pmh) & 
             !is.na(highest_level) & 
             !is.na(DemFH) & 
             !is.na(PDFH) & 
             !is.na(p26206) & 
             !is.na(p26260) & 
             !is.na(hyperlipidemia_pmh) & 
             !is.na(Hypertension_pmh) & 
             !is.na(obesity) & 
             !is.na(p22032_i0) & 
             !is.na(p31) & 
             !is.na(smoking) & 
             !is.na(p2020_i0))


##########Time before or during

cohort_of_interest_df <- subset(Overall, Overall$cohort_of_interest=="Yes")

Head_pre <- subset(cohort_of_interest_df, cohort_of_interest_df$Time_to_injury_1<cohort_of_interest_df$Age_at_recruitment)
Head_post <- subset(cohort_of_interest_df, cohort_of_interest_df$Time_to_injury_1>cohort_of_interest_df$Age_at_recruitment)

td_df_lag_before <- subset(td_df,  !td_df$eid%in%Head_post$eid)
Adj_2_td_df_lag_before <- subset(Adj_2_td_df_lag,  !Adj_2_td_df_lag$eid%in%Head_post$eid)

cox_fit_A <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = Adj_2_td_df_lag_before, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_A)

crude_incidence <- Adj_2_td_df_lag_before %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )

cox_fit_before <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury + Age_at_recruitment+alcohol+
                          +  Depression_pmh+p22189+Diabetes_pmh+highest_level+
                          +  DemFH+PDFH+p26206+p26260+
                          +  hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0
                        +  p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Adj_2_td_df_lag_before, id = eid, ties='efron' , timefix = TRUE)


summary(cox_fit_before)



Adj_2_td_df_lag_after <- subset(Adj_2_td_df_lag,  !Adj_2_td_df_lag$eid%in%Head_pre$eid)
td_df_lag_after <- subset(td_df,  !td_df$eid%in%Head_pre$eid)


cox_fit_A <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = td_df_lag_after, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_A)

crude_incidence <- td_df_lag_after %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )




cox_fit_after <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury + Age_at_recruitment+alcohol+
                         +  Depression_pmh+p22189+Diabetes_pmh+highest_level+
                         +  DemFH+PDFH+p26206+p26260+
                         +  hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0
                       +  p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Adj_2_td_df_lag_after, id = eid, ties='efron' , timefix = TRUE)


summary(cox_fit_after)

######Comm or GP only

#######################################################################################

Head_injured_Comm <- read.csv("~/Daniel_Whitehouse/Head_injured_cohort/output/complete_head_injured_Core_Comm_self.csv")
key <- Adj_2 %>% select(c(eid, Age_at_recruitment, p31, p2020_i0, visual_loss, hearing_loss, p22189, highest_level, p24003, p24006, p22032_i0, obesity, smoking,alcohol,DemFH, PDFH, hyperlipidemia_pmh,Hypertension_pmh, Diabetes_pmh, Cerebrovascular_pmh, Cardiovascular_pmh, Depression_pmh, p26206, p26260, start, stop, event,stop_dementia, event_dementia, time_in_study_dementia, stop_Parkinsons, event_Parkinsons, time_in_study_Parkinsons))
Key_neuro <- Outcome_dataset %>% select(eid, Time_to_Neurodegeneration, p40007_i0)

Head_injured_Comm <- merge(Head_injured_Comm, key, by="eid", all.y=T)
Head_injured_Comm <- merge(Head_injured_Comm, Key_neuro, by="eid", all.y=T)

Head_injured_Comm$head_injured_death <- ifelse(
  !is.na(Head_injured_Comm$Time_to_injury_1) & 
    !is.na(Head_injured_Comm$p40007_i0) & 
    abs(Head_injured_Comm$Time_to_injury_1 - Head_injured_Comm$p40007_i0) <= 0.25, 
  1,0
)


Head_injured_Comm$head_injury_start_1 <- ifelse(is.na(Head_injured_Comm$Time_to_injury_1), NA, Head_injured_Comm$Time_to_injury_1)

Head_injured_Comm$head_injury_start_1 <- ifelse(Head_injured_Comm$head_injured_death==1, NA, Head_injured_Comm$Time_to_injury_1)

Head_injured_Comm$head_injury_start_1 <- ifelse(Head_injured_Comm$head_injury_start_1 < Head_injured_Comm$Age_at_recruitment, 
                                                (Head_injured_Comm$Age_at_recruitment+0.01), 
                                                Head_injured_Comm$head_injury_start_1)

Head_injured_Comm$head_injury_start_1 <- Head_injured_Comm$head_injury_start_1 - Head_injured_Comm$Age_at_recruitment


Head_injured_Comm$Time_to_injury_1 <- as.numeric(Head_injured_Comm$Time_to_injury_1)
Head_injured_Comm$Time_to_Neurodegeneration <- as.numeric(Head_injured_Comm$Time_to_Neurodegeneration)

initial_na_counts <- sapply(1:1, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Head_injured_Comm[[time_col]]))
})

for (i in 1:1) {
  start_col <- paste("head_injury_start_", i, sep = "")
  time_col <- paste("Time_to_injury_", i, sep = "")
  condition <- Head_injured_Comm[[time_col]] > (Head_injured_Comm$Time_to_Neurodegeneration - 0.25)
  Head_injured_Comm[[start_col]][condition] <- NA
}

updated_na_counts <- sapply(1:1, function(i) {
  time_col <- paste("head_injury_start_", i, sep = "")
  sum(is.na(Head_injured_Comm[[time_col]]))
})


na_changes <- updated_na_counts - initial_na_counts
print(na_changes)

Simple_3_month_Comm <- Head_injured_Comm %>% select(c(eid, start, stop, event, head_injury_start_1, severity_1, head_injured_death))

Hosp_head <- subset(Simple_3_month_HES, !is.na(Simple_3_month_HES$head_injury_start_1))
Simple_3_month_Comm <- subset(Simple_3_month_Comm, !Simple_3_month_Comm$eid%in%Hosp_head$eid)

td_df_Comm <- tmerge(Simple_3_month_Comm, Simple_3_month_Comm, id = eid, outcome = event(stop, event))
td_df_Comm <- tmerge(td_df_Comm, Simple_3_month_Comm, id = eid, head_injury = tdc(head_injury_start_1))
td_df_Comm$head_injury <- as.factor(td_df_Comm$head_injury)
cox_fit_B_Comm <- coxph(formula = Surv(tstart, tstop, outcome) ~ head_injury, data = td_df_Comm, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_B_Comm)

export_cox_summary_to_csv(cox_fit_B_Comm, "cox_fit_GP_removed_Unadj_TBI.csv")

crude_incidence <- td_df_Comm %>%
  group_by(head_injury) %>%
  summarise(
    events = sum(outcome, na.rm = TRUE),                
    person_time = sum(tstop - tstart, na.rm = TRUE),    
    incidence_rate = events / person_time * 1000,       
    ci_lower = (events / person_time - 1.96 * sqrt(events) / person_time) * 1000, 
    ci_upper = (events / person_time + 1.96 * sqrt(events) / person_time) * 1000  
  )

#Adjustment set 1: Age_at_recruitment, p31, p2020_i0+visual_loss+hearing_loss+p24003+p24006, p22189, highest_level, p24003, p24006, p22032_i0, obesity, smoking, alcohol     

Adj_Comm <- merge(Adj_1_cov,Simple_3_month_Comm, by= "eid", all.y = T)
Adj_Comm <- merge(Adj_2_cov,Adj_Comm, by= "eid", all.y = T)

Adj_Comm_df_2 <- tmerge(Adj_Comm, Adj_Comm, id = eid, outcome = event(stop, event))
Adj_Comm_df_2 <- tmerge(Adj_Comm_df_2, Adj_Comm, id = eid, head_injury = tdc(head_injury_start_1))
Adj_Comm_df_2$head_injury <- as.factor(Adj_Comm_df_2$head_injury)

Adj_Comm_df_2$Age_group <- cut(Adj_Comm_df_2$Age_at_recruitment, breaks = c(-Inf, 50, 60, Inf), labels = c("<50", "50-60", ">60"))
cox_fit_Comm_effects_Comp <- coxph(formula = Surv(tstart, tstop, outcome) ~  head_injury + Age_at_recruitment+alcohol+
                                     Depression_pmh+p22189+Diabetes_pmh+highest_level+
                                     DemFH+PDFH+p26206+p26260+
                                     hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Adj_Comm_df_2, id = eid, ties='efron' , timefix = TRUE)
summary(cox_fit_Comm_effects_Comp)


filtered_data <- Adj_Comm %>%
  filter(  !is.na(visual_loss) &
             !is.na(hearing_loss) &
             !is.na(p24003) &
             !is.na(p24006) &
             !is.na(p22189) &
             !is.na(Age_at_recruitment) & 
             !is.na(alcohol) & 
             !is.na(Depression_pmh) & 
             !is.na(Diabetes_pmh) & 
             !is.na(highest_level) & 
             !is.na(DemFH) & 
             !is.na(PDFH) & 
             !is.na(p26206) & 
             !is.na(p26260) & 
             !is.na(hyperlipidemia_pmh) & 
             !is.na(Hypertension_pmh) & 
             !is.na(obesity) & 
             !is.na(p22032_i0) & 
             !is.na(p31) & 
             !is.na(smoking) & 
             !is.na(p2020_i0))

################Neurocognition

Neurocog_summ <- read.csv("~/Daniel_Whitehouse/Neurocognition/Output/Summary_UKB5_Neurocognition_Adjusted_longitudinal.csv")
Neurocog_summ <- subset(Neurocog_summ, Neurocog_summ$eid %in% Overall$eid)
Neurocog <- merge(Overall, Neurocog_summ, by = "eid", all.y = T)

Neurocog$head_injured_prior <- ifelse(Neurocog$Time_to_injury_1<Neurocog$Age_at_recruitment, "Yes", "No")
Neurocog$head_injured_prior <- ifelse(is.na(Neurocog$head_injured_prior), "No", Neurocog$head_injured_prior)
Neurocog$head_injured_prior <- as.factor(Neurocog$head_injured_prior)

Neurocog <- Neurocog %>%
  filter(!is.na(PC1_Score) |
           !is.na(UKB_RT_i0) | 
           !is.na(UKB_Pairs_Matching_i0) |
           !is.na(UKB_Prospective_Memory_i0) | 
           !is.na(UKB_Fluid_IQ_i0) |
           !is.na(UKB_Numeric_memory_i0))

Neurocog$head_injury_status <- with(Neurocog, ifelse(
  is.na(Time_to_injury_1) | Time_to_injury_1 > Age_at_recruitment, 
  "No head injury", 
  ifelse(
    Time_to_injury_1 < Age_at_recruitment & 
      (is.na(Time_to_injury_2) | Time_to_injury_2 > Age_at_recruitment),
    "Single head injury", 
    ifelse(
      Time_to_injury_2 < Age_at_recruitment,
      "2+ head injuries", 
      NA 
    )
  )
))

hist(Neurocog$PC1_Score)

####Association between head injury and summary neurocog score

model_one <- lm(PC1_Score ~ head_injured_prior, data = Neurocog)
summary(model_one)

model_two <- lm(PC1_Score ~ head_injured_prior + Age_at_recruitment+alcohol+
                  Depression_pmh+p22189+Diabetes_pmh+highest_level+
                  hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog)
summary(model_two)

summary_one <- summary(model_one)
coef_one <- summary_one$coefficients
head_injured_prior_one <- coef_one["head_injured_priorYes", ]


summary_two <- summary(model_two)
coef_two <- summary_two$coefficients
head_injured_prior_two <- coef_two["head_injured_priorYes", ]

summary_table <- data.frame(
  Model = c("Model 1", "Model 2"),
  Estimate = c(head_injured_prior_one["Estimate"], head_injured_prior_two["Estimate"]),
  Std_Error = c(head_injured_prior_one["Std. Error"], head_injured_prior_two["Std. Error"]),
  t_value = c(head_injured_prior_one["t value"], head_injured_prior_two["t value"]),
  p_value = c(head_injured_prior_one["Pr(>|t|)"], head_injured_prior_two["Pr(>|t|)"])
)

n_subjects_one <- nobs(model_one)
n_subjects_two <- nobs(model_two)


print(summary_table)
write.csv(summary_table, "summary_table_neurocog.csv", row.names = FALSE)

####multiple

Neurocog$head_injury_status <-as.factor(Neurocog$head_injury_status)
Neurocog$head_injury_status <- relevel(Neurocog$head_injury_status, ref = "No head injury")

# Fit the models
model_one <- lm(PC1_Score ~ head_injury_status, data = Neurocog)
model_data <- model.frame(model_one, na.action = na.omit)
table(model_data$head_injury_status)

model_two <- lm(PC1_Score ~ head_injury_status + Age_at_recruitment+alcohol+
                  Depression_pmh+p22189+Diabetes_pmh+highest_level+
                  hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog)
summary(model_two)
model_data <- model.frame(model_two, na.action = na.omit)
table(model_data$head_injury_status)




summary_one <- summary(model_one)
coef_one <- summary_one$coefficients


single_one <- coef_one["head_injury_statusSingle head injury", ]
multi_one <- coef_one["head_injury_status2+ head injuries", ]


summary_two <- summary(model_two)
coef_two <- summary_two$coefficients


single_two <- coef_two["head_injury_statusSingle head injury", ]
multi_two <- coef_two["head_injury_status2+ head injuries", ]


summary_table <- data.frame(
  Model = rep(c("Model 1", "Model 2"), each = 2),
  Injury_Status = rep(c("Single head injury", "2+ head injuries"), 2),
  Estimate = c(single_one["Estimate"], multi_one["Estimate"], 
               single_two["Estimate"], multi_two["Estimate"]),
  Std_Error = c(single_one["Std. Error"], multi_one["Std. Error"], 
                single_two["Std. Error"], multi_two["Std. Error"]),
  t_value = c(single_one["t value"], multi_one["t value"], 
              single_two["t value"], multi_two["t value"]),
  p_value = c(single_one["Pr(>|t|)"], multi_one["Pr(>|t|)"], 
              single_two["Pr(>|t|)"], multi_two["Pr(>|t|)"])
)


n_subjects_one <- nobs(model_one)
n_subjects_two <- nobs(model_two)

summary_table$N_Subjects <- rep(c(n_subjects_one, n_subjects_two), each = 2)


print(summary_table)
write.csv(summary_table, "summary_table_neurocog_multi.csv", row.names = FALSE)

print("Summary table has been saved as 'summary_table_neurocog_multi.csv'.")

model_two_data <- subset(Neurocog, !is.na(Neurocog$PC1_Score))

table(model_two_data$head_injury_status)

filtered_data <- Neurocog %>%
  filter(!is.na(PC1_Score) &
           !is.na(Age_at_recruitment) & 
           !is.na(alcohol) & 
           !is.na(Depression_pmh) & 
           !is.na(Diabetes_pmh) & 
           !is.na(highest_level) & 
           !is.na(hyperlipidemia_pmh) & 
           !is.na(Hypertension_pmh) & 
           !is.na(obesity) & 
           !is.na(p22032_i0) & 
           !is.na(p31) & 
           !is.na(smoking) & 
           !is.na(p2020_i0))

##Time between injury

Neurocog <- Neurocog %>%
  mutate(across(starts_with("Time_to_injury_"), 
                ~ ifelse(. > Age_at_recruitment, NA, .), 
                .names = "filtered_{col}"))

Neurocog$closest_injury <- Neurocog %>%
  select(starts_with("filtered_Time_to_injury_")) %>%
  apply(1, max, na.rm = TRUE)


Neurocog$closest_injury[is.infinite(Neurocog$closest_injury)] <- NA

Neurocog <- Neurocog %>%
  select(-starts_with("filtered_"))

Neurocog$head_prior_time <- Neurocog$Age_at_recruitment - Neurocog$closest_injury

Neurocog$injury_time_cat <- cut(Neurocog$head_prior_time, breaks = c(-Inf, 1,5,10, Inf),
                                labels = c("within one", "one to five", "five to ten", "ten or greater"),
                                included.lowest = T)

Neurocog$injury_time_cat <- if_else(is.na(Neurocog$injury_time_cat), "No head injury", Neurocog$injury_time_cat)

Neurocog$injury_time_cat <- factor(Neurocog$injury_time_cat, levels = c("No head injury", "within one", "one to five", "five to ten", "ten or greater"))

Linear_time <- lm(PC1_Score ~ injury_time_cat + Age_at_recruitment + p31 + highest_level, data = Neurocog)
summary(Linear_time)

summary_TIME <- summary(Linear_time)

coefficients_df <- as.data.frame(summary_TIME$coefficients)
coefficients_df <- coefficients_df[, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]

write.csv(coefficients_df, "Linear_neruocog_summary.csv", row.names = TRUE)


###Individual tests


##Reaction time - originally Mean time to correctly identify matches, so longer is bad
##UKB_Pairs_Matching - originallyNumber of incorrect matches in round, so bigger is bad
##UKB_Prospective_Memory - binary 1 is good
##UKB_Fluid_IQ - bigger is good
##UKB_Numeric_memory - max correct so big is good

##Reaction time

model_one_RT <- lm(UKB_RT_i0 ~ head_injured_prior, data = Neurocog)
summary(model_one_RT)

model_two_RT <- lm(UKB_RT_i0 ~ head_injured_prior + Age_at_recruitment+alcohol+
                     Depression_pmh+p22189+Diabetes_pmh+highest_level+
                     hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog)
summary(model_two_RT)

Neurocog$Predicted_RT <- predict(model_two, newdata = Neurocog)

##UKB_Pairs_Matching

model_one_PairM <- lm(UKB_Pairs_Matching_i0 ~ head_injured_prior, data = Neurocog)
summary(model_one_PairM)

model_two_PairM <- lm(UKB_Pairs_Matching_i0 ~ head_injured_prior + Age_at_recruitment+alcohol+
                        Depression_pmh+p22189+Diabetes_pmh+highest_level+
                        hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog)
summary(model_two_PairM)


##UKB_Prospective_Memory
Neurocog$UKB_Prospective_Memory_i0 <- as.factor(Neurocog$UKB_Prospective_Memory_i0)
model_one_ProsM_binary <- glm(UKB_Prospective_Memory_i0 ~ head_injured_prior, data = Neurocog, family = binomial)
summary(model_one_ProsM_binary)

install.packages("fmsb")
library(fmsb)

model_two_ProsM_binary <- glm(UKB_Prospective_Memory_i0 ~ head_injured_prior + Age_at_recruitment+alcohol+
                                Depression_pmh+p22189+Diabetes_pmh+highest_level+
                                hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, 
                              data = Neurocog, family = binomial)
summary(model_two_ProsM_binary)
NagelkerkeR2(model_two_ProsM_binary)
exp(coef(model_one_ProsM_binary))  
exp(coef(model_two_ProsM_binary)) 
exp(confint(model_two_ProsM_binary)) 
##UKB_Fluid_IQ

summary_table <- Neurocog %>%
  group_by(head_injured_prior) %>%
  summarise(
    mean_IQ = mean(UKB_Fluid_IQ_i0, na.rm = TRUE),
    sd_IQ = sd(UKB_Fluid_IQ_i0, na.rm = TRUE),
    count = n()
  )
summary_table
model_one_IQ <- lm(UKB_Fluid_IQ_i0 ~ head_injured_prior, data = Neurocog)
summary(model_one_IQ)

model_two_IQ <- lm(UKB_Fluid_IQ_i0 ~ head_injured_prior + Age_at_recruitment+alcohol+
                     Depression_pmh+p22189+Diabetes_pmh+highest_level+
                     hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog)
summary(model_two_IQ)


##UKB_Numeric_memory
Neurocog$UKB_Numeric_memory_i0 <- as.numeric(Neurocog$UKB_Numeric_memory_i0)
model_one_NM <- lm(UKB_Numeric_memory_i0 ~ head_injured_prior, data = Neurocog)
summary(model_one_NM)

model_two_NM <- lm(UKB_Numeric_memory_i0 ~ head_injured_prior + Age_at_recruitment+alcohol+
                     Depression_pmh+p22189+Diabetes_pmh+highest_level+
                     hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog)
summary(model_two_NM)


##UKB_5

model_one_UKB5 <- lm(PC1_Score ~ head_injured_prior, data = Neurocog)
summary(model_one_UKB5)

model_two_UKB5 <- lm(PC1_Score ~ head_injured_prior + Age_at_recruitment+alcohol+
                     Depression_pmh+p22189+Diabetes_pmh+highest_level+
                     hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog)
summary(model_two_UKB5)


extract_head_injury_results <- function(model, test_name, model_name) {

  results <- tidy(model) %>%
    filter(term == "head_injured_priorYes") %>%  
    select(term, estimate, std.error, statistic, p.value) %>%
    mutate(Test = test_name, Model = model_name)
  
  
  n_participants <- nobs(model) 
  results <- results %>%
    mutate(Participants = n_participants)
  
  return(results)
}

# Extract results for each model
results_list <- list(
  extract_head_injury_results(model_one_RT, "Reaction Time", "Model 1"),
  extract_head_injury_results(model_two_RT, "Reaction Time", "Model 2"),
  extract_head_injury_results(model_one_PairM, "Pairs Matching", "Model 1"),
  extract_head_injury_results(model_two_PairM, "Pairs Matching", "Model 2"),
  extract_head_injury_results(model_one_ProsM_binary, "Prospective Memory", "Model 1"),
  extract_head_injury_results(model_two_ProsM_binary, "Prospective Memory", "Model 2"),
  extract_head_injury_results(model_one_IQ, "Fluid IQ", "Model 1"),
  extract_head_injury_results(model_two_IQ, "Fluid IQ", "Model 2"),
  extract_head_injury_results(model_one_NM, "Numeric Memory", "Model 1"),
  extract_head_injury_results(model_two_NM, "Numeric Memory", "Model 2"),
  extract_head_injury_results(model_one_UKB5, "UKB 5", "Model 1"),
  extract_head_injury_results(model_two_UKB5, "UKB 5", "Model 2")
)

# Combine results into a single dataframe
final_results <- bind_rows(results_list) %>%
  select(Test, Model, estimate, std.error, statistic, p.value, Participants)  

final_results_mult <- subset(final_results, final_results$Model=="Model 2")
final_results_mult$Adj_P <- p.adjust(final_results_mult$p.value, method = "fdr")

# Save to CSV file
write.csv(final_results, "Head_Injury_Results_with_Participants.csv", row.names = FALSE)

# View results
View(final_results)


###Multiple


Neurocog$head_injury_status <- factor(Neurocog$head_injury_status, 
                                      levels = c("No head injury", "Single head injury", "2+ head injuries"))

##Reaction time

model_one_RT <- lm(UKB_RT_i0 ~ head_injury_status, data = Neurocog)
summary(model_one_RT)

model_two_RT <- lm(UKB_RT_i0 ~ head_injury_status + Age_at_recruitment+alcohol+
                     Depression_pmh+p22189+Diabetes_pmh+highest_level+
                     hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog)

summary(model_two_RT)

######UKB_Pairs_Matching

model_one_PairM <- lm(UKB_Pairs_Matching_i0 ~ head_injury_status, data = Neurocog)
summary(model_one_PairM)

model_two_PairM <- lm(UKB_Pairs_Matching_i0 ~ head_injury_status + Age_at_recruitment+alcohol+
                        Depression_pmh+p22189+Diabetes_pmh+highest_level+
                        hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog)
summary(model_two_PairM)
model_data <- model.frame(model_two_PairM, na.action = na.omit)
table(model_data$head_injury_status)
plot(model_two_PairM)

##UKB_Prospective_Memory

model_one_ProsM_binary <- glm(UKB_Prospective_Memory_i0 ~ head_injury_status, data = Neurocog, family = binomial)
summary(model_one_ProsM_binary)

model_two_ProsM_binary <- glm(UKB_Prospective_Memory_i0 ~ head_injury_status + Age_at_recruitment+alcohol+
                                Depression_pmh+p22189+Diabetes_pmh+highest_level+
                                hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, 
                              data = Neurocog, family = binomial)
summary(model_two_ProsM_binary)
NagelkerkeR2(model_two_ProsM_binary)

exp(coef(model_two_ProsM_binary))  
exp(confint(model_two_ProsM_binary)) 

##UKB_Fluid_IQ

model_one_IQ <- lm(UKB_Fluid_IQ_i0 ~ head_injury_status, data = Neurocog)
summary(model_one_IQ)

model_two_IQ <- lm(UKB_Fluid_IQ_i0 ~ head_injury_status + Age_at_recruitment+alcohol+
                     Depression_pmh+p22189+Diabetes_pmh+highest_level+
                     hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog)
summary(model_two_IQ)
model_data <- model.frame(model_two_IQ, na.action = na.omit)
table(model_data$head_injury_status)

##UKB_Numeric_memory
Neurocog$UKB_Numeric_memory_i0 <- as.numeric(Neurocog$UKB_Numeric_memory_i0)
model_one_NM <- lm(UKB_Numeric_memory_i0 ~ head_injury_status, data = Neurocog)
summary(model_one_NM)

model_two_NM <- lm(UKB_Numeric_memory_i0 ~ head_injury_status + Age_at_recruitment+alcohol+
                     Depression_pmh+p22189+Diabetes_pmh+highest_level+
                     hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog)
summary(model_two_NM)

model_data <- model.frame(model_two_NM, na.action = na.omit)
table(model_data$head_injury_status)

##UKB_5

model_one_UKB5 <- lm(PC1_Score ~ head_injury_status, data = Neurocog)
summary(model_one_UKB5)

model_two_UKB5 <- lm(PC1_Score ~ head_injury_status + Age_at_recruitment+alcohol+
                     Depression_pmh+p22189+Diabetes_pmh+highest_level+
                     hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog)
summary(model_two_UKB5)

model_data <- model.frame(model_two_UKB5, na.action = na.omit)
table(model_data$head_injury_status)

extract_head_injury_results_multi <- function(model, test_name, model_name) {
 
  results <- tidy(model) %>%
    filter(term %in% c("head_injury_status2+ head injuries", "head_injury_statusSingle head injury")) %>%
    select(term, estimate, std.error, statistic, p.value) %>%
    mutate(Test = test_name, Model = model_name)
  
  n_participants <- nobs(model)  
  results <- results %>%
    mutate(Participants = n_participants)
  
  return(results)
}

# Extract results for each model, including head injury status variables
results_list <- list(
  extract_head_injury_results_multi(model_one_RT, "Reaction Time", "Model 1"),
  extract_head_injury_results_multi(model_two_RT, "Reaction Time", "Model 2"),
  extract_head_injury_results_multi(model_one_PairM, "Pairs Matching", "Model 1"),
  extract_head_injury_results_multi(model_two_PairM, "Pairs Matching", "Model 2"),
  extract_head_injury_results_multi(model_one_ProsM_binary, "Prospective Memory", "Model 1"),
  extract_head_injury_results_multi(model_two_ProsM_binary, "Prospective Memory", "Model 2"),
  extract_head_injury_results_multi(model_one_IQ, "Fluid IQ", "Model 1"),
  extract_head_injury_results_multi(model_two_IQ, "Fluid IQ", "Model 2"),
  extract_head_injury_results_multi(model_one_NM, "Numeric Memory", "Model 1"),
  extract_head_injury_results_multi(model_two_NM, "Numeric Memory", "Model 2"),
  extract_head_injury_results_multi(model_one_UKB5, "UKB 5", "Model 1"),
  extract_head_injury_results_multi(model_two_UKB5, "UKB 5", "Model 2")
)

# Combine results into a single dataframe
final_results <- bind_rows(results_list) %>%
  select(Test, Model, term, estimate, std.error, statistic, p.value, Participants)  # Include 'term' column for head injury statuses

final_results_mult <- subset(final_results, final_results$Model=="Model 2")
final_results_mult$Adj_P <- p.adjust(final_results_mult$p.value, method = "fdr")

write.csv(final_results, "Individual_multiple_neurocog.csv", row.names = FALSE)

##linear test of trend

Neurocog$head_injury_status_numeric <- recode(Neurocog$head_injury_status, 
                                              "No head injury" = 0, 
                                              "1 head injury" = 1, 
                                              "2+ head injuries" = 2)

model_two_RT <- lm(UKB_RT_i0 ~ head_injury_status_numeric + Age_at_recruitment+alcohol+
                     Depression_pmh+p22189+Diabetes_pmh+highest_level+
                     hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog)
summary(model_two_RT)

model_two_PairM <- lm(UKB_Pairs_Matching_i0 ~ head_injury_status_numeric + Age_at_recruitment+alcohol+
                        Depression_pmh+p22189+Diabetes_pmh+highest_level+
                        hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog)
summary(model_two_PairM)

model_two_ProsM_binary <- glm(UKB_Prospective_Memory_i0 ~ head_injury_status_numeric + Age_at_recruitment+alcohol+
                                Depression_pmh+p22189+Diabetes_pmh+highest_level+
                                hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, 
                              data = Neurocog, family = binomial)
summary(model_two_ProsM_binary)

model_two_IQ <- lm(UKB_Fluid_IQ_i0 ~ head_injury_status_numeric + Age_at_recruitment+alcohol+
                     Depression_pmh+p22189+Diabetes_pmh+highest_level+
                     hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog)
summary(model_two_IQ)

model_two_NM <- lm(UKB_Numeric_memory_i0 ~ head_injury_status_numeric + Age_at_recruitment+alcohol+
                     Depression_pmh+p22189+Diabetes_pmh+highest_level+
                     hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog)
summary(model_two_NM)

model_two <- lm(PC1_Score ~ head_injury_status + Age_at_recruitment+alcohol+
                  Depression_pmh+p22189+Diabetes_pmh+highest_level+
                  hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog)
summary(model_two)

##Picture

install.packages("effects")
library(effects)

Neurocog_clean <- Neurocog[!is.na(Neurocog$PC1_Score) & 
                             !is.na(Neurocog$head_injury_status) & 
                             !is.na(Neurocog$Age_at_recruitment) & 
                             !is.na(Neurocog$alcohol) & 
                             !is.na(Neurocog$Depression_pmh) & 
                             !is.na(Neurocog$p22189) & 
                             !is.na(Neurocog$Diabetes_pmh) & 
                             !is.na(Neurocog$highest_level) & 
                             !is.na(Neurocog$hyperlipidemia_pmh) & 
                             !is.na(Neurocog$Hypertension_pmh) & 
                             !is.na(Neurocog$obesity) & 
                             !is.na(Neurocog$p22032_i0) & 
                             !is.na(Neurocog$p31) & 
                             !is.na(Neurocog$smoking) & 
                             !is.na(Neurocog$p2020_i0) & 
                             !is.na(Neurocog$visual_loss) & 
                             !is.na(Neurocog$hearing_loss) & 
                             !is.na(Neurocog$p24003) & 
                             !is.na(Neurocog$p24006), ]

model_two <- lm(PC1_Score ~ head_injury_status + Age_at_recruitment+alcohol+
                  Depression_pmh+p22189+Diabetes_pmh+highest_level+
                  hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog_clean)
summary(model_two)


predictions <- effect("head_injury_status", model_two)
predictions_df <- as.data.frame(predictions)
plot(predictions)

# Plot the adjusted predictions on the original scale

Adj_UKB <-ggplot(predictions_df, aes(x = head_injury_status, y = fit, color = head_injury_status, group = head_injury_status)) +
  geom_point(size = 3) +  
  geom_line(linewidth = 1, aes(group = 1), color="black", linetype = "dashed") +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.2, linewidth = 1) +  
  labs(title = "UKB-5", 
       y = "UKB 5 score",  
       x = "") +  
  theme_classic() +
  scale_color_manual(values = c("No head injury" = "#0072B2",    
                                "Single head injury" = "#009E73", 
                                "2+ head injuries" = "#D55E00")) +  
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") 
Adj_UKB



Neurocog_clean <- Neurocog[!is.na(Neurocog$UKB_RT_i0) & 
                             !is.na(Neurocog$head_injury_status) & 
                             !is.na(Neurocog$Age_at_recruitment) & 
                             !is.na(Neurocog$alcohol) & 
                             !is.na(Neurocog$Depression_pmh) & 
                             !is.na(Neurocog$p22189) & 
                             !is.na(Neurocog$Diabetes_pmh) & 
                             !is.na(Neurocog$highest_level) & 
                             !is.na(Neurocog$hyperlipidemia_pmh) & 
                             !is.na(Neurocog$Hypertension_pmh) & 
                             !is.na(Neurocog$obesity) & 
                             !is.na(Neurocog$p22032_i0) & 
                             !is.na(Neurocog$p31) & 
                             !is.na(Neurocog$smoking) & 
                             !is.na(Neurocog$p2020_i0) & 
                             !is.na(Neurocog$visual_loss) & 
                             !is.na(Neurocog$hearing_loss) & 
                             !is.na(Neurocog$p24003) & 
                             !is.na(Neurocog$p24006), ]

model_two_RT <- lm(UKB_RT_i0 ~ head_injury_status + Age_at_recruitment+alcohol+
                     Depression_pmh+p22189+Diabetes_pmh+highest_level+
                     hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog_clean)

summary(model_two_RT)

predictions <- effect("head_injury_status", model_two_RT)
predictions_df <- as.data.frame(predictions)


# Plot the adjusted predictions on the original scale

Adj_RT <- ggplot(predictions_df, aes(x = head_injury_status, y = fit, color = head_injury_status, group = head_injury_status)) +
  geom_point(size = 3) +  
  geom_line(linewidth = 1, aes(group = 1), color="black", linetype = "dashed") +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.2, linewidth = 1) +  
  labs(title = "Reaction Time", 
       y = "Reaction Time (log)",  
       x = "") +  
  theme_classic() +
  scale_color_manual(values = c("No head injury" = "#0072B2",    
                                "Single head injury" = "#009E73", 
                                "2+ head injuries" = "#D55E00")) +  
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
Adj_RT


##

Neurocog_clean <- Neurocog[!is.na(Neurocog$UKB_Pairs_Matching_i0) & 
                             !is.na(Neurocog$head_injury_status) & 
                             !is.na(Neurocog$Age_at_recruitment) & 
                             !is.na(Neurocog$alcohol) & 
                             !is.na(Neurocog$Depression_pmh) & 
                             !is.na(Neurocog$p22189) & 
                             !is.na(Neurocog$Diabetes_pmh) & 
                             !is.na(Neurocog$highest_level) & 
                             !is.na(Neurocog$hyperlipidemia_pmh) & 
                             !is.na(Neurocog$Hypertension_pmh) & 
                             !is.na(Neurocog$obesity) & 
                             !is.na(Neurocog$p22032_i0) & 
                             !is.na(Neurocog$p31) & 
                             !is.na(Neurocog$smoking) & 
                             !is.na(Neurocog$p2020_i0) & 
                             !is.na(Neurocog$visual_loss) & 
                             !is.na(Neurocog$hearing_loss) & 
                             !is.na(Neurocog$p24003) & 
                             !is.na(Neurocog$p24006), ]

model_two_PairM <- lm(UKB_Pairs_Matching_i0 ~ head_injury_status + Age_at_recruitment+alcohol+
                        Depression_pmh+p22189+Diabetes_pmh+highest_level+
                        hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog_clean)
summary(model_two_PairM)


predictions <- effect("head_injury_status", model_two_PairM)
predictions_df <- as.data.frame(predictions)
predictions_df$original_scale <- exp(predictions_df$fit)-1 
#predictions_df$original_scale <- signif(predictions_df$original_scale, digits = 1) # Reverse log transformation
predictions_df$se_original_scale <- exp(predictions_df$se)-1 
predictions_df$ci_lower <- predictions_df$original_scale - 1.96 * predictions_df$se_original_scale
predictions_df$ci_upper <- predictions_df$original_scale + 1.96 * predictions_df$se_original_scale


# Plot the adjusted predictions on the original scale

Adj_PM <-ggplot(predictions_df, aes(x = head_injury_status, y = fit, color = head_injury_status, group = head_injury_status)) +
  geom_point(size = 3) +  
  geom_line(linewidth = 1, aes(group = 1), color="black", linetype = "dashed") +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.2, linewidth = 1) +  
  labs(title = "Pairs Matching", 
       y = "Pairs Memory (log+1)",  
       x = "") +  
  theme_classic() +
  scale_color_manual(values = c("No head injury" = "#0072B2",    
                                "Single head injury" = "#009E73", 
                                "2+ head injuries" = "#D55E00")) +  
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
Adj_PM



odds_ratio_df <- data.frame(
  head_injury_status = c("No head injury", "Single head injury", "2+ head injuries"),
  OR = c(1, 0.8412483, 0.6205788),
  CI_lower = c(NA, 0.7480517, 0.4435570),
  CI_upper = c(NA, 0.9480392, 0.8802364),
  color = c("#0072B2", "#009E73", "#D55E00")   for single, Red for 2+
)

odds_ratio_df$head_injury_status <- factor(odds_ratio_df$head_injury_status, 
                                           levels = c("No head injury", "Single head injury", "2+ head injuries"))

# Create the plot
Adj_PrM <- ggplot(odds_ratio_df, aes(y = OR, x = head_injury_status, color = head_injury_status)) +
  geom_point(size = 3) +  # Dots for OR
  geom_line(linewidth = 1, aes(group = 1), color="black", linetype = "dashed") +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, linewidth =1) + 
  #geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +  # Red dashed reference line at OR = 1
  scale_color_manual(values = c("No head injury" = "#0072B2",
                                "Single head injury" = "#009E73",  
                                "2+ head injuries" = "#D55E00")) +  # Red
  labs(title = "Prospective Memory",
       x = "",
       y = "Odds Ratio of correct recall") +  
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") 



##

Neurocog_clean <- Neurocog[!is.na(Neurocog$UKB_Fluid_IQ_i0) & 
                             !is.na(Neurocog$head_injury_status) & 
                             !is.na(Neurocog$Age_at_recruitment) & 
                             !is.na(Neurocog$alcohol) & 
                             !is.na(Neurocog$Depression_pmh) & 
                             !is.na(Neurocog$p22189) & 
                             !is.na(Neurocog$Diabetes_pmh) & 
                             !is.na(Neurocog$highest_level) & 
                             !is.na(Neurocog$hyperlipidemia_pmh) & 
                             !is.na(Neurocog$Hypertension_pmh) & 
                             !is.na(Neurocog$obesity) & 
                             !is.na(Neurocog$p22032_i0) & 
                             !is.na(Neurocog$p31) & 
                             !is.na(Neurocog$smoking) & 
                             !is.na(Neurocog$p2020_i0) & 
                             !is.na(Neurocog$visual_loss) & 
                             !is.na(Neurocog$hearing_loss) & 
                             !is.na(Neurocog$p24003) & 
                             !is.na(Neurocog$p24006), ]

model_two_IQ <- lm(UKB_Fluid_IQ_i0 ~ head_injury_status + Age_at_recruitment+alcohol+
                     Depression_pmh+p22189+Diabetes_pmh+highest_level+
                     hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog_clean)



predictions <- effect("head_injury_status", model_two_IQ)
predictions_df <- as.data.frame(predictions)
predictions_df$ci_lower <- predictions_df$fit - 1.96 * predictions_df$se
predictions_df$ci_upper <- predictions_df$fit + 1.96 * predictions_df$se


# Plot the adjusted predictions on the original scale

Adj_IQ <-ggplot(predictions_df, aes(x = head_injury_status, y = fit, color = head_injury_status, group = head_injury_status)) +
  geom_point(size = 3) +  
  geom_line(linewidth = 1, aes(group = 1), color="black", linetype = "dashed") +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.2, linewidth = 1) +  
  labs(title = "Fluid Intelligence", 
       y = "Fluid IQ",  
       x = "") +  
  theme_classic() +
  scale_color_manual(values = c("No head injury" = "#0072B2",    
                                "Single head injury" = "#009E73", 
                                "2+ head injuries" = "#D55E00")) +  
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        legend.position = "none") 
Adj_IQ


#
Neurocog$UKB_Numeric_memory_i0 <- as.numeric(Neurocog$UKB_Numeric_memory_i0)
Neurocog_clean <- Neurocog[!is.na(Neurocog$UKB_Numeric_memory_i0) & 
                             !is.na(Neurocog$head_injury_status) & 
                             !is.na(Neurocog$Age_at_recruitment) & 
                             !is.na(Neurocog$alcohol) & 
                             !is.na(Neurocog$Depression_pmh) & 
                             !is.na(Neurocog$p22189) & 
                             !is.na(Neurocog$Diabetes_pmh) & 
                             !is.na(Neurocog$highest_level) & 
                             !is.na(Neurocog$hyperlipidemia_pmh) & 
                             !is.na(Neurocog$Hypertension_pmh) & 
                             !is.na(Neurocog$obesity) & 
                             !is.na(Neurocog$p22032_i0) & 
                             !is.na(Neurocog$p31) & 
                             !is.na(Neurocog$smoking) & 
                             !is.na(Neurocog$p2020_i0) & 
                             !is.na(Neurocog$visual_loss) & 
                             !is.na(Neurocog$hearing_loss) & 
                             !is.na(Neurocog$p24003) & 
                             !is.na(Neurocog$p24006), ]

model_two_NM <- lm(UKB_Numeric_memory_i0 ~ head_injury_status + Age_at_recruitment+alcohol+
                     Depression_pmh+p22189+Diabetes_pmh+highest_level+
                     hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = Neurocog_clean)
summary(model_two_NM)

predictions <- effect("head_injury_status", model_two_NM)
predictions_df <- as.data.frame(predictions)

# Plot the adjusted predictions on the original scale

Adj_NM <-ggplot(predictions_df, aes(x = head_injury_status, y = fit, color = head_injury_status, group = head_injury_status)) +
  geom_point(size = 3) +  
  geom_line(linewidth = 1, aes(group = 1), color="black", linetype = "dashed") +  
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.2, linewidth = 1) +  
  labs(title = "Numeric Memory", 
       y = "UKB Numeric Memory",  
       x = "") +  
  theme_classic() +
  scale_color_manual(values = c("No head injury" = "#0072B2",    
                                "Single head injury" = "#009E73", 
                                "2+ head injuries" = "#D55E00")) +  
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") 
Adj_NM


combined_plot <- plot_grid(Adj_RT, Adj_PM, Adj_PrM, Adj_IQ, Adj_NM, Adj_UKB, ncol = 2, 
                           rel_heights = c(1, 1, 1.2))
combined_plot
ggsave("Combined_Neurocog_Adj.tif", combined_plot, width = 7.5, height = 10)


#Longitudinal

##First time to test

Neurocog$p53_i1<-as.factor(Neurocog$p53_i1)
Neurocog$p53_i1<-strptime(Neurocog$p53_i1,format="%Y-%m-%d")
Neurocog$p53_i1<-as.Date(Neurocog$p53_i1,format="%Y-%m-%d")

Neurocog$p53_i2<-as.factor(Neurocog$p53_i2)
Neurocog$p53_i2<-strptime(Neurocog$p53_i2,format="%Y-%m-%d")
Neurocog$p53_i2<-as.Date(Neurocog$p53_i2,format="%Y-%m-%d")

Neurocog$p53_i3<-as.factor(Neurocog$p53_i3)
Neurocog$p53_i3<-strptime(Neurocog$p53_i3,format="%Y-%m-%d")
Neurocog$p53_i3<-as.Date(Neurocog$p53_i3,format="%Y-%m-%d")

Neurocog$Date_of_Birth_imputed<-as.factor(Neurocog$Date_of_Birth_imputed)
Neurocog$Date_of_Birth_imputed<-strptime(Neurocog$Date_of_Birth_imputed,format="%Y-%m-%d")
Neurocog$Date_of_Birth_imputed<-as.Date(Neurocog$Date_of_Birth_imputed,format="%Y-%m-%d")

Neurocog$Age_at_visit_1 <- difftime(Neurocog$p53_i1,Neurocog$Date_of_Birth_imputed,units=c("weeks"))
Neurocog$Age_at_visit_1 <-Neurocog$Age_at_visit_1/52.1429
Neurocog$Age_at_visit_1<-gsub("weeks","",as.character(Neurocog$Age_at_visit_1))
Neurocog$Age_at_visit_1<-as.numeric(Neurocog$Age_at_visit_1)

Neurocog$Age_at_visit_2 <- difftime(Neurocog$p53_i2,Neurocog$Date_of_Birth_imputed,units=c("weeks"))
Neurocog$Age_at_visit_2 <-Neurocog$Age_at_visit_2/52.1429
Neurocog$Age_at_visit_2<-gsub("weeks","",as.character(Neurocog$Age_at_visit_2))
Neurocog$Age_at_visit_2<-as.numeric(Neurocog$Age_at_visit_2)

Neurocog$Age_at_visit_3 <- difftime(Neurocog$p53_i3,Neurocog$Date_of_Birth_imputed,units=c("weeks"))
Neurocog$Age_at_visit_3 <-Neurocog$Age_at_visit_3/52.1429
Neurocog$Age_at_visit_3<-gsub("weeks","",as.character(Neurocog$Age_at_visit_3))
Neurocog$Age_at_visit_3<-as.numeric(Neurocog$Age_at_visit_3)

###


Neurocog$time_visit_0 <- 0
Neurocog$time_visit_1 <- Neurocog$Age_at_visit_1 - Neurocog$Age_at_recruitment
Neurocog$time_visit_2 <- Neurocog$Age_at_visit_2 - Neurocog$Age_at_recruitment
Neurocog$time_visit_3 <- Neurocog$Age_at_visit_3 - Neurocog$Age_at_recruitment


###Head during

Neurocog <- Neurocog %>%
  mutate(head_injury_during_study = case_when(
    Time_to_injury_1 >= Age_at_recruitment & 
      (Time_to_injury_1 <= Age_at_visit_2 | Time_to_injury_1 <= Age_at_visit_3) ~ 1,
    TRUE ~ 0
  ))

Neurocog <- Neurocog %>%
  mutate(Time_recruit_to_inj = ifelse(head_injury_during_study == 1, 
                                      Time_to_injury_1 - Age_at_recruitment, 
                                      NA))


install.packages("lme4")
library(lme4)
install.packages("lmerTest")
library(lmerTest)

##Reaction


df_long_RT <- Neurocog %>%
  filter(!is.na(UKB_RT_i1) | !is.na(UKB_RT_i2) | !is.na(UKB_RT_i3)) %>%
  select(eid, head_injured_prior, Time_to_injury_1, Age_at_recruitment, Age_at_visit_1, Age_at_visit_2, Age_at_visit_3, 
         UKB_RT_i1, UKB_RT_i2, UKB_RT_i3, UKB_RT_i0, alcohol,time_visit_0, time_visit_1, time_visit_2, time_visit_3,
         Depression_pmh, p22189, Diabetes_pmh, highest_level, head_injury_during_study, Time_to_injury_1,Time_recruit_to_inj,
         hyperlipidemia_pmh, Hypertension_pmh, obesity, p22032_i0, p31, smoking, p2020_i0, visual_loss, hearing_loss, p24003, p24006) %>%
  pivot_longer(cols = c(time_visit_0, time_visit_1,  time_visit_2, time_visit_3),
               names_to = "Timepoint", values_to = "Time") %>%
  mutate(UKB_RT = case_when(
    Timepoint == "time_visit_0" ~ UKB_RT_i0,
    Timepoint == "time_visit_1" ~ UKB_RT_i1,
    Timepoint == "time_visit_2" ~ UKB_RT_i2,
    Timepoint == "time_visit_3" ~ UKB_RT_i3
  ))

df_long_RT <- subset(df_long_RT, !is.na(df_long_RT$Time))
df_long_RT <- subset(df_long_RT, !is.na(df_long_RT$UKB_RT_i0))

df_long_RT <- df_long_RT[    !is.na(df_long_RT$UKB_RT) & 
                               !is.na(df_long_RT$head_injured_prior) & 
                               !is.na(df_long_RT$p22189) & 
                               !is.na(df_long_RT$highest_level) 
                             , ]

df_long_RT <- df_long_RT %>%
  mutate(Post_Injury = ifelse(is.na(Time_recruit_to_inj) | Time < Time_recruit_to_inj, 0, 1))



df_long_RT <- subset(df_long_RT, df_long_RT$head_injured_prior=="No")
df_long_RT <- subset(df_long_RT, !is.na(df_long_RT$Post_Injury))

###

df_long_RT$Age_centered <- as.numeric(scale(df_long_RT$Age_at_recruitment, center = TRUE, scale = TRUE))

lmm_reaction <- lmer(UKB_RT ~ Post_Injury*Time+highest_level+p31+Age_centered+(1|eid), data = df_long_RT)
summary(lmm_reaction)
print(lmm_reaction, correlation=T)


###UKB Pairs


df_long_Pairs <- Neurocog %>%
  filter(!is.na(UKB_Pairs_Matching_i1) | !is.na(UKB_Pairs_Matching_i2) | !is.na(UKB_Pairs_Matching_i3)) %>%
  select(eid, head_injured_prior, Time_to_injury_1, Age_at_recruitment, Age_at_visit_1, Age_at_visit_2, Age_at_visit_3, 
         UKB_Pairs_Matching_i1, UKB_Pairs_Matching_i2, UKB_Pairs_Matching_i3, UKB_Pairs_Matching_i0, alcohol,time_visit_0, time_visit_1, time_visit_2, time_visit_3,
         Depression_pmh, p22189, Diabetes_pmh, highest_level, head_injury_during_study, Time_to_injury_1,Time_recruit_to_inj,
         hyperlipidemia_pmh, Hypertension_pmh, obesity, p22032_i0, p31, smoking, p2020_i0, visual_loss, hearing_loss, p24003, p24006) %>%
  pivot_longer(cols = c(time_visit_0, time_visit_1,  time_visit_2, time_visit_3),
               names_to = "Timepoint", values_to = "Time") %>%
  mutate(UKB_Pairs = case_when(
    Timepoint == "time_visit_0" ~ UKB_Pairs_Matching_i0,
    Timepoint == "time_visit_1" ~ UKB_Pairs_Matching_i1,
    Timepoint == "time_visit_2" ~ UKB_Pairs_Matching_i2,
    Timepoint == "time_visit_3" ~ UKB_Pairs_Matching_i3
  ))

df_long_Pairs <- subset(df_long_Pairs, !is.na(df_long_Pairs$Time))
df_long_Pairs <- subset(df_long_Pairs, !is.na(df_long_Pairs$UKB_Pairs_Matching_i0))


df_long_Pairs <- df_long_Pairs[!is.na(df_long_Pairs$UKB_Pairs) & 
                                 !is.na(df_long_Pairs$head_injured_prior) & 
                                 !is.na(df_long_Pairs$p22189) & 
                                 !is.na(df_long_Pairs$highest_level) 
                               , ]

df_long_Pairs <- df_long_Pairs %>%
  mutate(Post_Injury = ifelse(is.na(Time_recruit_to_inj) | Time < Time_recruit_to_inj, 0, 1))


df_long_Pairs <- subset(df_long_Pairs, df_long_Pairs$head_injured_prior=="No")
df_long_Pairs <- subset(df_long_Pairs, !is.na(df_long_Pairs$Post_Injury))

###

df_long_Pairs$Age_centered <- as.numeric(scale(df_long_Pairs$Age_at_recruitment, center = TRUE, scale = TRUE))

lmm_pairs <- lmer(UKB_Pairs ~ Post_Injury*Time+highest_level+p31+Age_centered+(1|eid), data = df_long_Pairs)
summary(lmm_pairs)
print(lmm_pairs, correlation=T)


###UKB Prospective memory


df_long_ProsMem <- Neurocog %>%
  filter(!is.na(UKB_Prospective_Memory_i1) | !is.na(UKB_Prospective_Memory_i2) | !is.na(UKB_Prospective_Memory_i3)) %>%
  select(eid, head_injured_prior, Time_to_injury_1, Age_at_recruitment, Age_at_visit_1, Age_at_visit_2, Age_at_visit_3, 
         UKB_Prospective_Memory_i1, UKB_Prospective_Memory_i2, UKB_Prospective_Memory_i3, UKB_Prospective_Memory_i0, alcohol,time_visit_0, time_visit_1, time_visit_2, time_visit_3,
         Depression_pmh, p22189, Diabetes_pmh, highest_level, head_injury_during_study, Time_to_injury_1,Time_recruit_to_inj,
         hyperlipidemia_pmh, Hypertension_pmh, obesity, p22032_i0, p31, smoking, p2020_i0, visual_loss, hearing_loss, p24003, p24006) %>%
  pivot_longer(cols = c(time_visit_0, time_visit_1,  time_visit_2, time_visit_3),
               names_to = "Timepoint", values_to = "Time") %>%
  mutate(UKB_ProsMem = case_when(
    Timepoint == "time_visit_0" ~ UKB_Prospective_Memory_i0,
    Timepoint == "time_visit_1" ~ UKB_Prospective_Memory_i1,
    Timepoint == "time_visit_2" ~ UKB_Prospective_Memory_i2,
    Timepoint == "time_visit_3" ~ UKB_Prospective_Memory_i3
  ))

df_long_ProsMem <- subset(df_long_ProsMem, !is.na(df_long_ProsMem$Time))
df_long_ProsMem <- subset(df_long_ProsMem, !is.na(df_long_ProsMem$UKB_Prospective_Memory_i0))


df_long_ProsMem <- df_long_ProsMem[!is.na(df_long_ProsMem$UKB_ProsMem) & 
                                     !is.na(df_long_ProsMem$head_injured_prior) & 
                                     !is.na(df_long_ProsMem$p22189) & 
                                     !is.na(df_long_ProsMem$highest_level) 
                                   , ]

df_long_ProsMem <- df_long_ProsMem %>%
  mutate(Post_Injury = ifelse(is.na(Time_recruit_to_inj) | Time < Time_recruit_to_inj, 0, 1))


df_long_ProsMem <- subset(df_long_ProsMem, df_long_ProsMem$head_injured_prior=="No")
df_long_ProsMem <- subset(df_long_ProsMem, !is.na(df_long_ProsMem$Post_Injury))

###

df_long_ProsMem$Age_centered <- as.numeric(scale(df_long_ProsMem$Age_at_recruitment, center = TRUE, scale = TRUE))

lmm_pros <- lmer(UKB_ProsMem ~ Post_Injury*Time+highest_level+p31+Age_centered+(1|eid), data = df_long_ProsMem)
summary(lmm_pros)
print(lmm_pros, correlation=T)

##Fluid IQ


df_long_Fluid <- Neurocog %>%
  filter(!is.na(UKB_Fluid_IQ_i1) | !is.na(UKB_Fluid_IQ_i2) | !is.na(UKB_Fluid_IQ_i3)) %>%
  select(eid, head_injured_prior, Time_to_injury_1, Age_at_recruitment, Age_at_visit_1, Age_at_visit_2, Age_at_visit_3, 
         UKB_Fluid_IQ_i1, UKB_Fluid_IQ_i2, UKB_Fluid_IQ_i3, UKB_Fluid_IQ_i0, alcohol,time_visit_0, time_visit_1, time_visit_2, time_visit_3,
         Depression_pmh, p22189, Diabetes_pmh, highest_level, head_injury_during_study, Time_to_injury_1,Time_recruit_to_inj,
         hyperlipidemia_pmh, Hypertension_pmh, obesity, p22032_i0, p31, smoking, p2020_i0, visual_loss, hearing_loss, p24003, p24006) %>%
  pivot_longer(cols = c(time_visit_0, time_visit_1,  time_visit_2, time_visit_3),
               names_to = "Timepoint", values_to = "Time") %>%
  mutate(UKB_FluidIQ = case_when(
    Timepoint == "time_visit_0" ~ UKB_Fluid_IQ_i0,
    Timepoint == "time_visit_1" ~ UKB_Fluid_IQ_i1,
    Timepoint == "time_visit_2" ~ UKB_Fluid_IQ_i2,
    Timepoint == "time_visit_3" ~ UKB_Fluid_IQ_i3
  ))

df_long_Fluid <- subset(df_long_Fluid, !is.na(df_long_Fluid$Time))
df_long_Fluid <- subset(df_long_Fluid, !is.na(df_long_Fluid$UKB_Fluid_IQ_i0))


df_long_Fluid <- df_long_Fluid[!is.na(df_long_Fluid$UKB_Fluid_IQ_i0) & 
                                 !is.na(df_long_Fluid$head_injured_prior) & 
                                 !is.na(df_long_Fluid$p22189) & 
                                 !is.na(df_long_Fluid$highest_level) 
                               , ]

df_long_Fluid <- df_long_Fluid %>%
  mutate(Post_Injury = ifelse(is.na(Time_recruit_to_inj) | Time < Time_recruit_to_inj, 0, 1))


df_long_Fluid <- subset(df_long_Fluid, df_long_Fluid$head_injured_prior=="No")
df_long_Fluid <- subset(df_long_Fluid, !is.na(df_long_Fluid$Post_Injury))

###

df_long_Fluid$Age_centered <- as.numeric(scale(df_long_Fluid$Age_at_recruitment, center = TRUE, scale = TRUE))

lmm_fluid <- lmer(UKB_Fluid_IQ_i0 ~ Post_Injury*Time+highest_level+p31+Age_centered+(1|eid), data = df_long_Fluid)
summary(lmm_fluid)
print(lmm_fluid, correlation=T)

residuals <- resid(lmm_fluid)
fitted_values <- fitted(lmm_fluid)

plot(fitted_values, residuals, 
     main = "Residuals vs. Fitted Values", 
     xlab = "Fitted Values", 
     ylab = "Residuals")
abline(h = 0, col = "red")









##Cohort with summary neurocog prior to first injury

Neurocog_summ <- subset(Neurocog, !is.na(Neurocog$UKB_RT_i0))
Neurocog_summ <- subset(Neurocog_summ, !is.na(Neurocog_summ$UKB_RT_i1) | !is.na(Neurocog_summ$UKB_RT_i2)| !is.na(Neurocog_summ$UKB_RT_i3))
Neurocog_summ_no_head <- subset(Neurocog_summ, Neurocog_summ$head_injured_prior=="No")

Neurocog_summ_no_head$head_injury_between <- ifelse(
  is.na(Neurocog_summ_no_head$Time_to_injury_1), 0,  # If Time_to_injury_1 is NA, assign 0
  ifelse(
    Neurocog_summ_no_head$Time_to_injury_1 >= Neurocog_summ_no_head$Age_at_recruitment & 
      (Neurocog_summ_no_head$Time_to_injury_1 <= Neurocog_summ_no_head$Age_at_visit_1 | 
         Neurocog_summ_no_head$Time_to_injury_1 <= Neurocog_summ_no_head$Age_at_visit_2 |
         Neurocog_summ_no_head$Time_to_injury_1 <= Neurocog_summ_no_head$Age_at_visit_3), 
    "Yes", "No"
  )
)

Neurocog_summ_no_head$head_injury_between <- as.factor(Neurocog_summ_no_head$head_injury_between)

#Neurocog_summ_no_head <- subset(Neurocog_summ_no_head, Neurocog_summ_no_head$head_injury_between=="Yes")
df_long <- Neurocog_summ_no_head %>%
  select(eid, head_injury_between, Time_to_injury_1, Age_at_recruitment, Age_at_visit_1, Age_at_visit_2, Age_at_visit_3, 
         UKB_RT_i1, UKB_RT_i2, UKB_RT_i3, UKB_RT_i0) %>%
  pivot_longer(cols = c(Age_at_recruitment,Age_at_visit_1,  Age_at_visit_2, Age_at_visit_3),
               names_to = "Timepoint", values_to = "Age") %>%
  mutate(UKB_RT = case_when(
    Timepoint == "Age_at_recruitment" ~ UKB_RT_i0,
    Timepoint == "Age_at_visit_1" ~ UKB_RT_i1,
    Timepoint == "Age_at_visit_2" ~ UKB_RT_i2,
    Timepoint == "Age_at_visit_3" ~ UKB_RT_i3,
  ))

df_long <- subset(df_long, !is.na(df_long$Age))
df_long <- subset(df_long, !is.na(df_long$UKB_RT))

df_long <- df_long %>%
  arrange(eid, Age) %>%
  group_by(eid) %>%
  mutate(head_injury_status = ifelse(!is.na(Time_to_injury_1) & Time_to_injury_1 <= Age, 1, 0))

install.packages("lme4")
library(lme4)

lmm_reaction <- lmer(UKB_RT ~ Age * head_injury_status + (1 | eid), data = df_long)
summary(lmm_reaction)

df_long$Predicted_RT <- predict(lmm_reaction)

ggplot(df_long, aes(x = Age, color = as.factor(head_injury_status))) +
  geom_smooth(aes(y = UKB_RT), method = "loess", se = FALSE, linetype = "solid", size = 1) +  # Actual
  geom_smooth(aes(y = Predicted_RT), method = "loess", se = FALSE, linetype = "dashed", size = 1) +  # Predicted
  labs(title = "Actual vs. Predicted Reaction Time",
       x = "Age",
       y = "Reaction Time (ms)",
       color = "Head Injury Status") +
  theme_minimal()


########################################################With imputation##########################################################################################


###Adjusted for:
#A) Age_at_recruitment, p31, p2020_i0, p22189, highest_level, p24003, p24006, p22032_i0, obesity, p20117_i0,p20116_i0
#B) Age_at_recruitment, p31, p2020_i0, p22189, highest_level, p24003, p24006, p22032_i0, obesity, p20117_i0,p20116_i0, Hypertension_pmh, Diabetes_pmh, Stroke_pmh, Cardiovascular_pmh, depression_pmh, p26206, p26260, p21000_i0
For_imputation <- Overall
cov <- For_imputation %>% select(c(eid, Age_at_recruitment,alcohol,
                                     Depression_pmh,p22189,Diabetes_pmh,highest_level,
                                     DemFH,PDFH,p26206,p26260,
                                     hyperlipidemia_pmh,Hypertension_pmh,obesity,p22032_i0,p31,smoking,p2020_i0,visual_loss,hearing_loss,
                                     p24003,p24006))
cov <- subset(cov, !is.na(cov$Age_at_recruitment))

Outcome_for_imputation <- Simple_3_months_lag %>% select(c(c("eid",
                                                             "start",             
                                                             "stop" ,             
                                                             "event",              
                                                             "head_injury_start_1", 
                                                             "head_injury_start_2",
                                                             "time_in_study")))

Imp <- merge(cov, Outcome_for_imputation, by= "eid", all.x = T)

Imp$Age_at_recruitment <- as.numeric(Imp$Age_at_recruitment)
Imp$alcohol <- factor(Imp$alcohol, levels = c("Never", "Previous", "Light", "Moderate", "Heavy"))
Imp$Depression_pmh <- as.factor(Imp$Depression_pmh)
Imp$p22189 <- as.numeric(Imp$p22189)
Imp$Diabetes_pmh <- as.factor(Imp$Diabetes_pmh)
Imp$highest_level <- factor(Imp$highest_level, levels = c("Primary","Secondary", "Training","College"))
Imp$DemFH <- as.factor(Imp$DemFH)
Imp$PDFH <- as.factor(Imp$PDFH)
Imp$p26206 <- as.numeric(Imp$p26206)
Imp$p26260 <- as.numeric(Imp$p26260)
Imp$hyperlipidemia_pmh <- as.factor(Imp$hyperlipidemia_pmh)
Imp$Hypertension_pmh <- as.factor(Imp$Hypertension_pmh)
Imp$obesity <- as.factor(Imp$obesity)
Imp$p22032_i0 <- factor(Imp$p22032_i0, levels = c("high","moderate", "low"))
Imp$p31 <- as.factor(Imp$p31)
Imp$smoking <- factor(Imp$smoking, levels = c("Never", "Previous", "Occasional", "Frequent"))
Imp$p2020_i0 <- as.factor(Imp$p2020_i0)
Imp$visual_loss <- as.factor(Imp$visual_loss)
Imp$hearing_loss <- as.factor(Imp$hearing_loss)
Imp$p24003 <- as.numeric(Imp$p24003)
Imp$p24006 <- as.numeric(Imp$p24006)

#Imputation using MICE 

init = mice(Imp, maxit=0) 
meth = init$method
predM = init$predictorMatrix

predM[, colnames(predM) %in% c("eid",
                               "start",             
                               "stop" ,             
                               "event",              
                               "head_injury_start_1", 
                               "head_injury_start_2",
                               "time_in_study")] <- 0

meth[names(meth) %in% c( "eid",
                         "start",             
                         "stop" ,             
                         "event",              
                         "head_injury_start_1",  
                         "head_injury_start_2",
                         "time_in_study")] <- ""

meth[names(meth) %in% c( "Age_at_recruitment",
                         "p22189",
                         "p26206",
                         "p26260", 
                         "p24003",
                         "p24006" )] <- "pmm"

meth[names(meth) %in% c( "Depression_pmh",      
                         "Diabetes_pmh",
                         "DemFH",
                         "PDFH",
                         "hyperlipidemia_pmh",
                         "Hypertension_pmh",
                         "obesity",
                         "p31",
                         "p2020_i0",
                         "visual_loss",
                         "hearing_loss"
                          )] <- "logreg"

meth[names(meth) %in% c( "highest_level",
                         "p22032_i0",
                         "smoking",
                         "alcohol")] <- "polyreg"

imputed_Data <- mice(Imp, m=5, predictorMatrix=predM, seed=500, method = meth)

imputed_data_long <- complete(imputed_Data, action = 'long', include = TRUE)

write.csv(imputed_data_long, "imputed_Data_TBI.csv", row.names = FALSE)

#imputed_data_long <- read.csv("~/Daniel_Whitehouse/Neurodegeneration/output/TBI/Imputed/imputed_Data_TBI.csv", header = T)

imputed_Data <- as.mids(imputed_data_long)


###Primary analysis#############################################################################################################

complete_data_list <- lapply(1:5, function(i) complete(imputed_Data, i))

###Full adjustment Imputed

common_ids <- Simple_3_months_lag$eid

# Subset imputed dataframes to remove those not in 1 year#########################################################

complete_data_list_1_year <- lapply(1:5, function(i) {
  imputed_data_i <- complete(imputed_Data, i)
  imputed_data_subset <- imputed_data_i[imputed_data_i$eid %in% common_ids, ]
  return(imputed_data_subset)
})

tmerge_list <- lapply(complete_data_list_1_year, function(df) {
  df_tmerge <- tmerge(df, df, id = eid, outcome = event(stop, event))
  df_tmerge <- tmerge(df_tmerge, df, id = eid, head_injury = tdc(head_injury_start_1))
  df_tmerge$head_injury <- as.factor(df_tmerge$head_injury)
  return(df_tmerge)
})

cox_models_1year_Adj_1 <- lapply(tmerge_list, function(df_tmerge) {
  coxph(Surv(tstart, tstop, outcome) ~ head_injury + Age_at_recruitment+alcohol+
          Depression_pmh+p22189+Diabetes_pmh+highest_level+
          DemFH+PDFH+p26206+p26260+
          hyperlipidemia_pmh+Hypertension_pmh+obesity+p22032_i0+p31+smoking+p2020_i0+visual_loss+hearing_loss+p24003+p24006, data = df_tmerge, ties='efron' , timefix = FALSE)
})

summary(pool(cox_models_1year_Adj_1))
pooled_results <- pool(cox_models_1year_Adj_1)
summary_results <- summary(pooled_results)
summary_results$HR <- exp(summary_results$estimate)
summary_results$lower_CI <- exp(summary_results$estimate - 1.96 * summary_results$std.error)
summary_results$upper_CI <- exp(summary_results$estimate + 1.96 * summary_results$std.error)
summary_results$p_value <- summary_results$p.value

imputed_result_total <- as.data.frame(summary_results[, c("term", "HR", "lower_CI", "upper_CI", "p_value")])
write.csv(imputed_result_total, "Imputed_demographics_adjusted_1_year_lag_TBI.csv")

cox_models_1year_Adj_2 <- lapply(tmerge_list, function(df_tmerge) {
  coxph(Surv(tstart, tstop, outcome)  ~ head_injury + strata(Age_group) + p31 + p2020_i0 + p22189 + highest_level + p24003 + p24006 + p22032_i0 + obesity+ smoking + alcohol + Hypertension_pmh + Diabetes_pmh + Cerebrovascular_pmh + Cardiovascular_pmh + Depression_pmh +p26206+ p26260, data = df_tmerge, ties='efron' , timefix = FALSE)
})

summary(pool(cox_models_1year_Adj_2))
pooled_results <- pool(cox_models_1year_Adj_2)
summary_results <- summary(pooled_results)
summary_results$HR <- exp(summary_results$estimate)
summary_results$lower_CI <- exp(summary_results$estimate - 1.96 * summary_results$std.error)
summary_results$upper_CI <- exp(summary_results$estimate + 1.96 * summary_results$std.error)
summary_results$p_value <- summary_results$p.value

imputed_result_total <- as.data.frame(summary_results[, c("term", "HR", "lower_CI", "upper_CI", "p_value")])
write.csv(imputed_result_total, "Imputed_total_adjustment_3_month_lag_TBI.csv")

