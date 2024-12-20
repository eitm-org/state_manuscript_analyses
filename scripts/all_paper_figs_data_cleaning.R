#BEFORE YOU CAN RUN THIS YOU WILL NEED:
#1) to set up a rclone remote to dropbox called "dropbox"
#to see filepaths you may have to change, search "^"

# Package Installations & Set-Up --------------------------------------------------------
if(!require(tidyverse)){
  install.packages("tidyverse")
  library(tidyverse)
}

if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}

if(!require(ggbeeswarm)){
  install.packages("ggbeeswarm")
  library(ggbeeswarm)
}

if(!require(patchwork)){
  install.packages("patchwork")
  library(patchwork)
}

if(!require(viridis)){
  install.packages("viridis")
  library(viridis)
}

if(!require(cowplot)){
  install.packages("cowplot")
  library(cowplot)
}

if(!require(stringr)){
  install.packages("stringr")
  library(stringr)
}

if(!require(janitor)){
  install.packages("janitor")
  library(janitor)
}

if(!require(gridExtra)){
  install.packages("gridExtra")
  library(gridExtra)
}

if(!require(stringr)){
  install.packages("stringr")
  library(stringr)
}

if(!require(sjPlot)){
  install.packages("sjPlot")
  library(sjPlot)
}

if(!require(sjtable2df)){
  install.packages("sjtable2df")
  library(sjtable2df)
}

if(!require(lme4)){
  install.packages("lme4")
  library(lme4)
}

if(!require(readxl)){
  install.packages("readxl")
  library(readxl)
}

if(!require(purrr)){
  install.packages("purrr")
  library(purrr)
}

# Working directory should be repo root folder (^)
new_wd <- file.path("~", "project", "state_manuscript_analyses")
setwd(new_wd)

source(file.path("scripts", "all_paper_figs_functions.R"))

######################################################
#
#        DATA PATHS
#
######################################################
# Contrived samples
if(!require(reticulate)){
  install.packages("reticulate")
  library(reticulate)
}
source_python('post_processing/scripts/constants.py', envir = parent.frame())
contrived <- contrived_path

# sbs cosmic fit data
SBS_fp <- file.path("input_data/fig3/SBS_false_positives.csv")
SBS_filtered <- file.path("input_data/fig3/SBS_post_filtering.csv")
SBS_subject <- file.path("input_data/fig3/SBS_pre_filtering.csv")
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
######################################################
#
#        FIG 1A: CONTRIVED SAMPLES DATA CLEANING
#
######################################################

samples <- c("EIBS-002GA_19226_1", "EIBS-002GB_19227_1", "EIBS-002GC_19228_1",
             "EIBS-002GD_19229_1", "EIBS-002GA_19226_2", "EIBS-002GB_19227_2",
             "EIBS-002GC_19228_2", "EIBS-002GD_19229_2", "EIBS-002GB_19227_3",
             "EIBS-002GC_19228_3", "EIBS-002GD_19229_4", "EIBS-002GB_19227_4",
             "EIBS-002GC_19228_4", "EIBS-002GD_19229_3")

# Search for files that contain specified sample IDs
files0 <- list.files(contrived, full.names = TRUE)
files0 <- files0[sapply(files0, function(files) some(samples, function(x) grepl(x, files)))]
files0 <- files0[!grepl(".idx", files0)]

# Make input data folder to store VCFs
if (!dir.exists(file.path("input_data", "fig1", "recovery_pct"))){
  dir.create(file.path("input_data", "fig1", "recovery_pct"))
}

samp100files <- list.files(contrived, full.names = TRUE, recursive = TRUE)
samp100files <- samp100files[grepl("_germline.ann.germ.af.filtered3.vcf.gz", samp100files)]

for (i in 1:length(samp100files)) {
  system(paste("cp", samp100files[[i]], file.path("input_data", "fig1", "recovery_pct")))
  s100filename <- basename(samp100files[[i]])
  system(paste("bcftools index -f", file.path("input_data", "fig1", "recovery_pct", s100filename)))
  samp100files[[i]] <- file.path("input_data", "fig1", "recovery_pct", s100filename)
}

# Get intersection, compress, and copy ground truth vcf
system(paste("bcftools isec -p", file.path("input_data", "fig1", "recovery_pct", "s100isec"), samp100files[[1]], samp100files[[2]]))
gt_vcf <- file.path("input_data", "fig1", "recovery_pct", "s100isec", "0002.vcf")
system(paste("bgzip -f", gt_vcf))
gt_vcf <- paste(gt_vcf, ".gz", sep = "")
system(paste("cp ", gt_vcf, "input_data/fig1/recovery_pct/s100isec/0002_copy.vcf.gz"))
gt_vcf2 <- file.path("input_data/fig1/recovery_pct/s100isec/0002_copy.vcf.gz")
system(paste("bcftools index -f", gt_vcf))

results <- data.frame(sample = as.character(),
                      filename = as.character(),
                      sample_type = as.character(),
                      len0 = as.numeric(),
                      len1 = as.numeric(),
                      len2 = as.numeric(),
                      len3 = as.numeric())

files0 <- files0[!grepl("EIBS-002GA_19226", files0)]
files <- c(files0, samp100files, gt_vcf2)

# Iterate through files to get intersection with bt474.vcf
for (file in files) {
  
  if (!grepl("EIBS-002GA_19226", file) & file != gt_vcf2) {
    system(paste("cp", file, file.path("input_data", "fig1", "recovery_pct")))
    filename <- basename(file)
    isecdir <- paste(strsplit(filename, "\\.")[[1]][[1]], "vs100", sep = "_")
    # Zip the file using bcftools
    system(paste("bgzip -f", file.path("input_data", "fig1", "recovery_pct", filename)))
    zipfile <- file.path("input_data", "fig1", "recovery_pct", paste(filename, ".gz", sep = ""))
  } else if (file == gt_vcf2) {
    filename <- basename(file)
    isecdir <- paste(strsplit(filename, "\\.")[[1]][[1]], "vs100", sep = "_")
    zipfile <- file
  } else {
    filename <- basename(file)
    zipfile <- file
  }
  
  contrived_map <- c('EIBS-002GA' = "100% sample", 'EIBS-002GB' = "1% sample", 'EIBS-002GC' = "5% sample", 'EIBS-002GD' = "10% sample")
  
  system(paste("bcftools index -f", zipfile))
  # Getting complement of this file using bt474.vcf "ground truth" file
  system(paste("bcftools isec -p", file.path("input_data", "fig1", "recovery_pct", isecdir), zipfile, gt_vcf))
  sff_strt <- regexpr("EIBS-[[:digit:]]{3}[[:alpha:]]{2}_[[:digit:]]{5}_[[:digit:]]", filename)[[1]]
  sff_end <- sff_strt + attr(regexpr("EIBS-[[:digit:]]{3}[[:alpha:]]{2}_[[:digit:]]{5}_[[:digit:]]", filename), "match.length") - 1
  sampfromfile <- substr(filename, sff_strt, sff_end)
  if (sampfromfile != "") {
    sample_type <- contrived_map[[substring(sampfromfile, 0, 10)]]
  }
  if (file == gt_vcf2) {
    sample_type <- "100% intersection"
  }

  vcfs <- list.files(file.path("input_data", "fig1", "recovery_pct", isecdir))
  vcfs <- sort(vcfs[grepl(".vcf", vcfs)])
  vcflens <- c()
  for (i in vcfs) {
    vcf_file <- file.path("input_data", "fig1", "recovery_pct", isecdir, i)
    system("touch vcflens.txt")
    system(paste("bcftools view -H", vcf_file, "| wc -l >> vcflens.txt"))
    vcflens[[length(vcflens) + 1]] <- as.numeric(read.delim("vcflens.txt", header = FALSE)[[1]])
    system("rm vcflens.txt")
  }
  
  row <- c(sampfromfile, filename, sample_type, vcflens)
  names(row) <- names(results)
  results <- rbind(results, row)
  
}

# f1 Score Calculation
results <- results %>%
  mutate(tp = len2,
         #0002 = true positive
         #0000 = false positive
         fp = len0,
         fn = len1) %>%
  mutate(precision = tp / (tp + fp),
         error_rate = 1 - precision,
         recall = tp / (tp + fn)*100,
         f1 = 2*(precision*recall) / (precision + recall),
         sample_purity = as.numeric(str_remove(str_extract(sample_type, "[[:digit:]]+%"), "%")))

write_rds(results, file.path("input_data", "fig1", "recovery_pct", "isec_results_vs100sample.rds"))

########################################################
#
#     SUPPLEMENTAL FIG: SBS MUTSIGS FOR SUBJECTS VS FP
#
########################################################

fp_df <- read.csv(SBS_fp)
post_df <- read.csv(SBS_filtered)
pre_df <- read.csv(SBS_subject)

# make sampled dataframes
snv_ratios <- bind_rows(post_df, fp_df)
other_cols <- find_sbs_others(snv_ratios, 0.05, 0.85)

set.seed(1) # random seed
sample_df <- bind_rows(
  post_df[sample(nrow(post_df), size = 100), ], 
  fp_df[sample(nrow(fp_df), size = 20), ]
)

# significant vs others
sig_mutsigs <- subset(sample_df, select = c('Data.Type', 'Sample.Names', other_cols))
other_mutsigs <- sample_df[ , !(names(sample_df) %in% other_cols)]

sig_mutsigs$other <- rowSums(other_mutsigs[-1:-2], na.rm=TRUE)
other_mutsigs$other <- rowSums(sig_mutsigs[3:(length(sig_mutsigs) - 1)], na.rm=TRUE)

# normalize for each sample
sig_mutsigs[-1:-2] <- sig_mutsigs[-1:-2] / rowSums(sig_mutsigs[-1:-2], na.rm=TRUE)
other_mutsigs[-1:-2] <- other_mutsigs[-1:-2] / rowSums(other_mutsigs[-1:-2], na.rm=TRUE)

sig_mutsigs <- sig_mutsigs[-1] %>% pivot_longer(-'Sample.Names')
other_mutsigs <- other_mutsigs[-1] %>% pivot_longer(-'Sample.Names') 

# save data
write_csv(sig_mutsigs, file.path("input_data", "fig3", "sig_mutsigs.csv"))
write_csv(other_mutsigs, file.path("input_data", "fig3", "other_mutsigs.csv"))

############################################################
#
#        WAFFLE / CATEGORICAL PLOT OF SBS MUTSIGS
#
############################################################

cat_mutsigs <- bind_rows(fp_df, post_df, pre_df)
cat_mutsigs <- cat_mutsigs %>% 
  pivot_longer(cols = starts_with("SBS"),
    names_to = "signature",
    values_to = "ratio",
    values_drop_na = TRUE)
  
cat_mutsigs$ratio <- cat_mutsigs$ratio / 100

write_csv(cat_mutsigs, file.path("input_data", "fig3", "cat_mutsigs.csv"))

############################################################
#
#        FIG 5: TMB PLOT
#
############################################################

activities_df <- read_csv(file.path("input_data", "fig5", "activities.csv"))
plotting_data <- read_csv(file.path("input_data", "fig5", "data_for_tmb.csv"))

activities_df <- activities_df %>% group_by(Types) %>% arrange(Mut_burden, .by_group=TRUE)
activities_df$index <- rep(c(1:length(activities_df$Types[activities_df$Types == 'SBS1'])), 40)
activities_df <- activities_df[activities_df$Mut_burden != 0, ]

x_labs <- list()
for (name in plotting_data$names) {
  x_labs <- append(x_labs, median(activities_df[activities_df$Types == name, ]$index))
}
plotting_data$x_labs <- as.numeric(x_labs)

write_csv(activities_df, file.path("input_data", "fig5", "long_activities.csv"))
write_csv(plotting_data, file.path("input_data", "fig5", "long_data_for_tmb.csv"))

############################################################
#
#        FIG 6: LQMM REGRESSION PLOTS
#
############################################################

# LQMM AGE SUMMARY
lqmm_age_summaries <- lqmm_for_age(new_wd)
lqmm_age_coeff <- filter(lqmm_age_summaries, Type == 'draw_age' & Variable != 'distance_to_transcript_count_adj3')
lqmm_age_coeff <- lqmm_age_coeff %>% rename('lqmm.p' = 'Pr(>|t|)', Std.Error = 'Std. Error')

lqmm_age_coeff <- lqmm_age_coeff[!is.na(lqmm_age_coeff[['lqmm.p']]), ] # drop null values

# save
if (!dir.exists(file.path("input_data", "fig6"))) {
  dir.create(file.path("input_data", "fig6"))
}
write_csv(lqmm_age_coeff, './input_data/fig6/lqmm_age_data.csv')

## ANOVA SUMMARY
anova_summary_age <- anova_for_age(new_wd)
anova_summary_age <- anova_summary_age %>% rename('anova.p' = 'p-value')
anova_summary_age <- anova_summary_age[!is.na(anova_summary_age[['anova.p']]), ] # drop NA p-values

write_csv(anova_summary_age, file.path("input_data", "fig6", "anova_summary_age.csv"))

## MAKE TABLE
# separate intercept and slope data
intercepts <- lqmm_age_summaries[lqmm_age_summaries$Type == '(Intercept)', ]
slopes <- lqmm_age_summaries[lqmm_age_summaries$Type == 'draw_age', ]

make_table(intercepts, slopes, new_wd)




