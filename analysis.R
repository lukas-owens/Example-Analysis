############################################################
# Project: Breast Cancer Survival Analysis
# Author: Lukas Owens, Fred Hutchinson Cancer Center
# Date: 7/11/2024
# Description: This example project analyzes a cohort of women
#    treated for localized breast cancer with surgery. Various
#    baseline variables are provided, as is all-cause survival.
#    The analysis produces Kaplan-Meier survival curves by 
#    receptor status, and runs a Cox proportional hazards 
#    regression of survival on the available covariates.
#
############################################################

library(tidyverse)
library(survival)
library(survminer)
library(gtsummary)

datestamp <- '2024-07-11'
data_path <- 'data/breast_cancer_raw_2024-07-11.csv'

############################################################
# Import and clean data set for use in analysis. The variables included are:
#     - id: Patient ID
#     - year_diag: Year of diagnosis with breast cancer
#     - dob: Date of birth
#     - tumor_size: Size of tumor (in mm)
#     - nodes: Number of lymph nodes with cancer
#     - receptors: ER/PGR receptor status
#     - time: Time of death or censoring
#     - status: Death indicator
############################################################
process_data <- function(path) {
    dset <- read_csv(path, skip = 1, col_names = c('id',
                                                   'year_diag',
                                                   'dob',
                                                   'tumor_size',
                                                   'grade',
                                                   'nodes',
                                                   'receptors',
                                                   'time',
                                                   'status'),
                     col_types = 'ciDccicdi')
    
    # Calculate age at diagnosis (assuming diagnoses at mid-year)
    dset <- dset |> mutate(age_diag = ymd(paste0(year_diag, '0701')) - dob)
    dset <- dset |> mutate(age_diag = time_length(age_diag, 'years'))
    dset <- dset |> select(-dob)
    
    # Turn "unknown" values into explicit NAs
    dset <- dset |> mutate(grade = ifelse(grade == 'Unknown', NA, grade))
    
    # Convert categorical data to factors
    dset <- dset |> mutate(tumor_size = factor(tumor_size, levels = c('<=20', '20-50', '>50')))
    dset <- dset |> mutate(grade = factor(grade, levels = c('2', '3')))
    
    # Separate receptors into two columns
    dset <- dset |> mutate(er = ifelse(receptors %in% c('ER positive', 'Double positive'), 1, 0))
    dset <- dset |> mutate(pgr = ifelse(receptors %in% c('PGR positive', 'Double positive'), 1, 0))
    dset <- dset |> mutate(across(c(er, pgr), as.integer))
    
    # Convert time from days into years
    dset <- dset |> mutate(time = time / 365.25)
    
    dset
}

############################################################
# Summarize the baseline characteristics of the variables of
# interest for this study
############################################################
summarize_baseline_chars <- function(data, saveit = FALSE) {
    tbl <- data |> select(age_diag, tumor_size, grade, nodes, er, pgr)  
    tbl <- tbl_summary(tbl, label = list('age_diag' = 'Age at Diagnosis',
                                         'tumor_size' = 'Tumor size (mm)',
                                         'grade' = 'Grade',
                                         'nodes' = 'Number of nodes',
                                         'er' = 'Estrogen receptor positive',
                                         'pgr' = 'Progesterone receptor positive'))
    if(saveit) {
        tbl |> as_flex_table() |> 
            flextable::save_as_docx(path = str_glue('output/table_baseline_features_{datestamp}.docx'))
    }
    tbl
}

############################################################
# Fit a Cox proportional hazards model for time to death 
# with hormone status and other prognostic factors
############################################################
run_survival_regression <- function(data, saveit = FALSE) {
    formula <- Surv(time, status) ~ er + pgr + age_diag + tumor_size + grade + nodes
    mod <- coxph(formula, data)
    if(saveit) {
        tbl <- mod |> tbl_regression(label = list('age_diag' = 'Age at Diagnosis',
                                                  'tumor_size' = 'Tumor size (mm)',
                                                  'grade' = 'Grade',
                                                  'nodes' = 'Number of nodes',
                                                  'er' = 'Estrogen receptor positive',
                                                  'pgr' = 'Progesterone receptor positive'),
                                     exponentiate = TRUE) 
        tbl |> as_flex_table() |> 
            flextable::save_as_docx(path = str_glue('output/table_regression_{datestamp}.docx'))
    }
    mod
}

############################################################
# Calculate survival curves by receptor stats for ER and PGR
############################################################
plot_survival_curves <- function(data, receptor, saveit = FALSE) {
    formula <- as.formula(paste('Surv(time, status) ~', receptor))
    km <- surv_fit(formula, data = data)
    g <- ggsurvplot(km, data = data, censor = FALSE, risk.table = TRUE)
    if(saveit) {
        png(str_glue('output/figure_survival_er_{datestamp}.png'), 
            width = 500, height = 550, res = 96)
        print(g)
        dev.off()
    }
    g
}


# Run analysis
breast_cancer_data <- process_data(data_path)
save(breast_cancer_data, file = 'data/breast_cancer_processed_2024-07-11.RData')
summarize_baseline_chars(breast_cancer_data, saveit = TRUE)
run_survival_regression(breast_cancer_data, saveit = TRUE)
plot_survival_curves(breast_cancer_data, 'er', saveit = TRUE)
plot_survival_curves(breast_cancer_data, 'pgr', saveit = TRUE)
