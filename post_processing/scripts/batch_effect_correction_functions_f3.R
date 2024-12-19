## Import Library ---------

# remotes::install_github("stopsack/batchtma") # if batchtma not installed
library(remotes)
library(batchtma)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(DBI)
library(RSQLite)
library(Gmisc)
library(dplyr)

if(!require(reticulate)){
  install.packages("reticulate")
  library(reticulate)
}
setwd('~/projects/state_manuscript_analyses/post_processing')
source_python('scripts/constants.py')

global = read.csv("aggregate_data/aggregate_snv_patient_global_f3.csv")

## Merge and Preprocess -------
data = read.csv(sample_prep_table_path)
data <- data %>%
  rename(
    EIBS = Sample
  )

data = data[data$Sample_type == 'subject',]
data = data %>% filter(if_all(Median_Coverage, ~ . >= 20))
data$LSK = str_split_fixed(data$Ligation_Sequencing_Kit_Batch., " ", 2)[,1]
data = distinct(data)
merged_qc_global <- merge(global, data, by = "EIBS") 

## Adjust All Columns -------
#' Adjust batch effect for all columns
#' 
#' @param data Data for batch effect correction
#' @param keyword keyword of the column variable for batch effect correction
#' @param batch Batch for batch effect correction
#' @param confounder Confounder for batch correction
#' @param saveFile Choose whether to save file to current directory
#' @returns A new csv file with batch effect adjusted
#' @examples adjust_all_batch(merged_qc_rm, count, Date, q75_QSCORE_40HR, TRUE)
#' merged_qc_global_adj<-adjust_all_batch(merged_qc_global, count, Date, q75_QSCORE_40HR, FALSE)
adjust_all_batch <- function(data, keyword, batch, confounder, saveFile) {
  split_data <- unlist(strsplit(deparse(substitute(data)), split = "_"))
  df_name <- split_data[length(split_data)]
  data_adj <- data[, names(get(df_name))]
  
  # Convert keyword to a string
  keyword <- deparse(substitute(keyword))
  for (col in names(data)) {
    if (grepl(keyword, col) && !grepl("_min_", col) && !grepl("_max_", col) && +
        !grepl("read_count", col) && !grepl("passed", col) && !startsWith(col, "count_"))
    {
      
      temp <- try({
        data %>%
          adjust_batch(markers = col, batch = {{ batch }}, method = "standardize", confounders = {{ confounder }})
      }, silent = TRUE)
      
      # Check for errors
      if (class(temp) == "try-error") {
        print(paste("Error in adjusting batch for column:", col))
        next
      }
      
      adjusted_column <- setdiff(names(temp), names(data))
      
      # only keep the final result
      if (length(adjusted_column) > 1) {
        adjusted_column <- tail(adjusted_column, n = 1)
      }
      
      # merge into the adjusted dataframe
      if(col %in% names(data_adj)) {
        temp <- temp %>% dplyr::select(!!sym(adjusted_column))
      } else {
        temp <- temp %>% dplyr::select(!!sym(adjusted_column), !!sym(col))
      }
      data_adj <- cbind(data_adj, temp)
    }
  }
  
  if (saveFile){
    filename <- paste0(deparse(substitute(data)), "_",deparse(substitute(confounder)),"_adj.csv")
    write.csv(data_adj, filename, row.names = FALSE) 
  }
  
  return(data_adj)
}

# Adjust Global
merged_qc_global_adj<-adjust_all_batch(merged_qc_global, count, LSK, Median_Coverage, F)
write.csv(merged_qc_global_adj, "aggregate_data/global_batch_adj_f3.csv", row.names = FALSE) 

## Plot Data ---------

#' Plot raw and adjusted data side by side
#' 
#' @param data Data for plotting
#' @param marker Variable need to be batch adjusted
#' @param batch Batch for batch effect correction
#' @param confounder Confounder for batch correction
#' @returns A side to side plot (left: unadjusted; right: adjusted)
#' @example plot_side_to_side(merged_qc_rm, snv_count, Date, q75_QSCORE_40HR)
plot_side_to_side <- function(data, marker, batch, confounder) {
  marker_str <- deparse(substitute(marker))
  adj_marker <- sym(paste0(marker_str, "_adj4"))
  corrected <- data %>%
    adjust_batch(markers = {{ marker }}, batch = {{ batch }}, 
                 method = "ipw", confounders = {{ confounder }}) %>%
    plot_batch(marker = {{ adj_marker }}, batch = {{ batch }}, color = {{ confounder }})
  uncorrected <- data %>% plot_batch(marker = {{ marker }}, batch = {{ batch }}, color = {{ confounder }})
  side_to_side <- uncorrected + corrected
  return(side_to_side)
}
