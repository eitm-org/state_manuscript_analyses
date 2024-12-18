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

USER = Sys.getenv(c("USER"))

## Set Work Directory --------
setwd(pathJoin("/home", USER, "/projects/STATE_QC/data"))

chronqc_db_dir = pathJoin('/home', USER, 'dropbox/chronqc_db')
system(paste0('rclone copy dropbox:"EITM-Oxford Nanopore Team/STATE experiments/Derived Data Analysis/ChronQC/chronqc_db" ', chronqc_db_dir))
## Load Data ---------
con <- dbConnect(SQLite(), "/home/xchen@okta-oci.eitm.org/dropbox/chronqc_db/chronqc.stats.sqlite")
as.data.frame(dbListTables(con))
data <- dbReadTable(con, 'chronqc_stats_data')
data$active_CHANNEL_40HR = data$pore_CHANNEL_40HR + data$sequencing_CHANNEL_40HR
data$LSK = str_split_fixed(data$Ligation_Sequencing_Kit_Batch., " ", 2)[,1]
data = distinct(data)

global = read.csv("/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_snv_patient_global_f1.csv")
chrom = read.csv("/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_snv_patient_chrom_f1.csv")
mbps = read.csv("/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_snv_patient_Mbps_f1.csv")

## Merge and Preprocess -------

# merge dataset and exclude run 002 - 004
data <- data %>%
  rename(
    EIBS = Sample
  )

ids_wanted = c('PAU25537',
               'PAU13632',
               'PAU24514',
               'PAU24660',
               'PAU25467',
               'PAU24089',
               'PAU25134',
               'PAU25155',
               'PAU23576',
               'PAU13428',
               'PAU25526',
               'PAU25178',
               'PAU24905',
               'PAU14686',
               'PAU24824',
               'PAU24938',
               'PAU24406',
               'PAU24062',
               'PAU24986')

data <- data %>%
  rename(
    EIBS = Sample
  )
data = data %>% filter(if_all(Median_Coverage, ~ . >= 20)) %>%
  filter(if_all(EIBS, ~ . != "EIBS-001TC_02816_3_1_1_20231031" & . != 'EIBS-001VB_03588_9_1_1_20221021'))  %>%
  filter(! flowcell_id %in% ids_wanted)


merged_qc_global <- merge(global, data, by = "EIBS") 

merged_qc_chrom <- merge(chrom, data, by = "EIBS")
merged_qc_mbps <- merge(mbps, data, by = "EIBS") 
## Split Mbps and Chrom Dataframe
# by the 'chr' column
list_qc_chrom <- split(merged_qc_chrom, merged_qc_chrom$chr)

#list of mbps by 'chr', 'start' and 'end'
merged_qc_mbps$group <- paste(merged_qc_mbps$chr, merged_qc_mbps$start, merged_qc_mbps$end, sep="_")
list_qc_mbps <- split(merged_qc_mbps, merged_qc_mbps$group)

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
      print(col)
      
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
      print(adjusted_column)
      
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

# 10 Hr
# Adjust Global
merged_qc_global_10Hr_adj<-adjust_all_batch(merged_qc_global, count, LSK, Median_Coverage, F)
write.csv(merged_qc_global_10Hr_adj, "/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/global_batch_adj_f1.csv", row.names = FALSE) 

# Adjust Chrom 
adjusted_dfs_chrom <- lapply(list_qc_chrom, function(df_chrom) {
  adjust_all_batch(df_chrom, count, LSK, Median_Coverage, FALSE)
}) # (lapply function name must contain "_chrom_")

all_columns <- unique(unlist(lapply(adjusted_dfs_chrom, names)))

adjusted_dfs_chrom <- lapply(adjusted_dfs_chrom, function(df) {
  missing_cols <- setdiff(all_columns, names(df))
  df[missing_cols] <- NA
  return(df)
})
merged_qc_chrom_10Hr_adj <- do.call(rbind, adjusted_dfs_chrom)
write.csv(merged_qc_chrom_10Hr_adj, "/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/chrom_batch_adj_f1.csv", row.names = FALSE) 

# Adjust mbps
adjusted_dfs_mbps <- lapply(list_qc_mbps, function(df_mbps) {
  adjust_all_batch(df_mbps, count, LSK, Median_Coverage, FALSE)
}) # (lapply function name must contain "_mbps_")

all_columns <- unique(unlist(lapply(adjusted_dfs_mbps, names)))

adjusted_dfs_mbps <- lapply(adjusted_dfs_mbps, function(df) {
  missing_cols <- setdiff(all_columns, names(df))
  df[missing_cols] <- NA
  return(df)
})
merged_qc_mbps_10Hr_adj <- do.call(rbind, adjusted_dfs_mbps)
write.csv(merged_qc_mbps_10Hr_adj, "/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/mbps_batch_adj_f1.csv", row.names = FALSE) 


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
