#BEFORE YOU CAN RUN THIS YOU WILL NEED:
#1) to set up a rclone remote to dropbox called "dropbox"

#to see filepaths you may have to change, search "^"

library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
library(patchwork)
#library(kableExtra)
library(viridis)
library(cowplot)
library(stringr)
library(janitor)
library(gridExtra)
library(sjPlot)
library(sjtable2df)
library(lme4)
library(readxl)
library(purrr)

new_wd <- file.path("~", "project", "STATE_analyses")
setwd(new_wd)

source(file.path("R", "all_paper_figs_functions.R"))

######################################################
#
#        DATA PATHS
#
######################################################
#contrived samples from xchen directory
xchen_contrived <- file.path("/data/scratch/xchen/STATE_vcfs_f3_region_filtered/full_genome")

#contrived samples from vaidhy's directory
#these are the 100% contrived samples-- vaidhy had to reconfigure the filtering for these, which is why they are in his directory
vaidhy_contrived <- file.path("/data/scratch/vmahaganapathy/contrived")

#xchen data for contrived snv metrics plot
xchen_contrived_snv_metrics <- file.path("/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/contrived_snv_metrics.csv")

# sbs cosmic fit data
SBS_fp <- file.path("input_data/fig3/SBS_false_positives.csv")
SBS_filtered <- file.path("input_data/fig3/SBS_post_filtering.csv")
SBS_subject <- file.path("input_data/fig3/SBS_pre_filtering.csv")

#fig5 snvs
fig5_snvs <- file.path("/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/input_data/fig5_joined_snv_clinical.csv")

#demographics data for csv
demo_redcap <- file.path("Redcap Cloud STATE", "STATE_Enrollment.csv")
cohort_data <- file.path("STATE Draw Event Data", "STATE Draws - Deidentified_cohorts.xlsx")
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

#search for filenames that contain these sample numbers in /data/scratch/xchen/STATE_vcfs_f3_region_filtered/full_genome directory
files0 <- list.files(xchen_contrived, full.names = TRUE)
#sapply() loops through filenames in the directory to...
#apply some(), which tells you if each individual filename has one of the contrived sample sample numbers in it
files0 <- files0[sapply(files0, function(files) some(samples, function(x) grepl(x, files)))]
files0 <- files0[!grepl(".idx", files0)]

#make input data folder to store VCFs
if (!dir.exists(file.path("input_data", "fig1", "recovery_pct"))){
  dir.create(file.path("input_data", "fig1", "recovery_pct"))
}

#set the sample 100 file as
samp100files <- list.files(vaidhy_contrived, full.names = TRUE, recursive = TRUE)
samp100files <- samp100files[grepl("_germline.ann.germ.af.filtered3.vcf.gz", samp100files)]
#copy it over and index
for (i in 1:length(samp100files)) {
  system(paste("cp", samp100files[[i]], file.path("input_data", "fig1", "recovery_pct")))
  s100filename <- basename(samp100files[[i]])
  system(paste("bcftools index -f", file.path("input_data", "fig1", "recovery_pct", s100filename)))
  samp100files[[i]] <- file.path("input_data", "fig1", "recovery_pct", s100filename)
}

#now get their intersection
system(paste("bcftools isec -p", file.path("input_data", "fig1", "recovery_pct", "s100isec"), samp100files[[1]], samp100files[[2]]))
#and use that as the ground truth
gt_vcf <- file.path("input_data", "fig1", "recovery_pct", "s100isec", "0002.vcf")
#compress
system(paste("bgzip -f", gt_vcf))
gt_vcf <- paste(gt_vcf, ".gz", sep = "")
#copy gt_vcf over to another one
system(paste("cp ", gt_vcf, "input_data/fig1/recovery_pct/s100isec/0002_copy.vcf.gz"))
gt_vcf2 <- file.path("input_data/fig1/recovery_pct/s100isec/0002_copy.vcf.gz")
#get index
system(paste("bcftools index -f", gt_vcf))

#make a dataframe to store your results
results <- data.frame(sample = as.character(),
                      filename = as.character(),
                      sample_type = as.character(),
                      len0 = as.numeric(),
                      len1 = as.numeric(),
                      len2 = as.numeric(),
                      len3 = as.numeric())
#samp100 files are filtered differently
#remove 100% purity samples from the main file directory
#and add the ones that were filtered correctly
files0 <- files0[!grepl("EIBS-002GA_19226", files0)]
files <- c(files0, samp100files, gt_vcf2)

#loop through files and get intersections with bt474.vcf
for (file in files) {
  
  #if the file is from one of the sample files, you won't need to copy, zip, etc it
  if (!grepl("EIBS-002GA_19226", file) & file != gt_vcf2) {
    system(paste("cp", file, file.path("input_data", "fig1", "recovery_pct")))
    filename <- basename(file)
    isecdir <- paste(strsplit(filename, "\\.")[[1]][[1]], "vs100", sep = "_")
    #first, zip the file the way bcftools wants you to
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
  #change the name of the sample column
  #next, index the file
  system(paste("bcftools index -f", zipfile))
  # #get the complement of this file with the bt474.vcf "ground truth" file
  system(paste("bcftools isec -p", file.path("input_data", "fig1", "recovery_pct", isecdir), zipfile, gt_vcf))
  sff_strt <- regexpr("EIBS-[[:digit:]]{3}[[:alpha:]]{2}_[[:digit:]]{5}_[[:digit:]]", filename)[[1]]
  sff_end <- sff_strt + attr(regexpr("EIBS-[[:digit:]]{3}[[:alpha:]]{2}_[[:digit:]]{5}_[[:digit:]]", filename), "match.length") - 1
  sampfromfile <- substr(filename, sff_strt, sff_end)
  sample_type <- samples_df[samples_df$`Fragmented DNA aliquot ID` == sampfromfile, "Contrived sample type"][[1]]
  if (file == gt_vcf2) {
    sample_type <- "100% intersection"
  }
  #loop through the vcf files in the isecdir you just made
  vcfs <- list.files(file.path("input_data", "fig1", "recovery_pct", isecdir))
  vcfs <- sort(vcfs[grepl(".vcf", vcfs)])
  vcflens <- c()
  for (i in vcfs) {
    # print(i)
    vcf_file <- file.path("input_data", "fig1", "recovery_pct", isecdir, i)
    #use touch to throw results of bcftools view -H into a text file
    system("touch vcflens.txt")
    system(paste("bcftools view -H", vcf_file, "| wc -l >> vcflens.txt"))
    #read text file into r
    vcflens[[length(vcflens) + 1]] <- as.numeric(read.delim("vcflens.txt", header = FALSE)[[1]])
    system("rm vcflens.txt")
  }
  
  #throw those results, + filename, sample_type and sampfromfile into a dataframe row
  row <- c(sampfromfile, filename, sample_type, vcflens)
  names(row) <- names(results)
  results <- rbind(results, row)
  
}

#calculate f1 score to plot!
results <- results %>%
  #0002 = true positive (intersection between both)
  mutate(tp = len2,
         #0000 = false positive (all the calls that we called in our sample that weren’t in ground truth)
         # our sequence minus bt474
         fp = len0,
         fn = len1) %>%
  mutate(precision = tp / (tp + fp),
         error_rate = 1 - precision,
         recall = tp / (tp + fn)*100,
         f1 = 2*(precision*recall) / (precision + recall),
         sample_purity = as.numeric(str_remove(str_extract(sample_type, "[[:digit:]]+%"), "%")))

write_rds(results, file.path("input_data", "fig1", "recovery_pct", "isec_results_vs100sample.rds"))

######################################################
#
#        CONTRIVED SNV METRICS
#
######################################################

if (!dir.exists(file.path("input_data", "fig1", "recovery_figs"))){
  dir.create(file.path("input_data", "fig1", "recovery_figs"))
}
system(paste("cp '", xchen_contrived_snv_metrics, "' '", file.path("input_data", "fig1", "recovery_figs", "contrived_snv_metrics.csv"), "'", sep = ""))


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
# write_csv(sig_mutsigs, file.path("input_data", "fig3", "sig_mutsigs.csv"))
# write_csv(other_mutsigs, file.path("input_data", "fig3", "other_mutsigs.csv"))

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
#        FIG 4: TMB PLOT
#
############################################################

activities_df <- read_csv(file.path("input_data", "fig4", "activities.csv"))
plotting_data <- read_csv(file.path("input_data", "fig4", "data_for_tmb.csv"))

activities_df <- activities_df %>% group_by(Types) %>% arrange(Mut_burden, .by_group=TRUE)
activities_df$index <- rep(c(1:length(activities_df$Types[activities_df$Types == 'SBS1'])), 40)
activities_df <- activities_df[activities_df$Mut_burden != 0, ]

x_labs <- list()
for (name in plotting_data$names) {
  x_labs <- append(x_labs, median(activities_df[activities_df$Types == name, ]$index))
}
plotting_data$x_labs <- as.numeric(x_labs)

write_csv(activities_df, file.path("input_data", "fig4", "activities.csv"))
write_csv(plotting_data, file.path("input_data", "fig4", "data_for_tmb.csv"))

############################################################
#
#        FIG 5: LQMM REGRESSION PLOTS
#
############################################################

# LQMM AGE SUMMARY
lqmm_age_summaries <- lqmm_for_age(new_wd)
lqmm_age_coeff <- filter(lqmm_age_summaries, Type == 'draw_age' & Variable != 'distance_to_transcript_count_adj3')
lqmm_age_coeff <- lqmm_age_coeff %>% rename('lqmm.p' = 'Pr(>|t|)', Std.Error = 'Std. Error')

lqmm_age_coeff <- lqmm_age_coeff[!is.na(lqmm_age_coeff[['lqmm.p']]), ] # drop null values

# save
if (!dir.exists(file.path("input_data", "fig5"))) {
  dir.create(file.path("input_data", "fig5"))
}
write_csv(lqmm_age_coeff, './input_data/fig5/lqmm_age_data.csv')

## ANOVA SUMMARY
anova_summary_age <- anova_for_age(new_wd)
anova_summary_age <- anova_summary_age %>% rename('anova.p' = 'p-value')
anova_summary_age <- anova_summary_age[!is.na(anova_summary_age[['anova.p']]), ] # drop NA p-values

write_csv(anova_summary_age, file.path("input_data", "fig5", "anova_summary_age.csv"))

## MAKE TABLE
# separate intercept and slope data
intercepts <- lqmm_age_summaries[lqmm_age_summaries$Type == '(Intercept)', ]
slopes <- lqmm_age_summaries[lqmm_age_summaries$Type == 'draw_age', ]

make_table(intercepts, slopes, new_wd)

############################################################
#
#        demographics csv
#
############################################################
system(paste("rclone copy 'dropbox:", demo_redcap, "' '", file.path("input_data", "demo"), "'", sep = ""))
system(paste("rclone copy 'dropbox:", cohort_data, "' '", file.path("input_data", "cohorts"), "'", sep = ""))

ppids <- c('9HZC16EL0', '9HZC17EM0', '9HZC17EMO', '9HZC17LGG', '9HZC18LFG',
           '9HZC18LEG', '9HZC18EMI', '9HZC19EI0', '9HZC19EL0', '9HZC19EJ0',
           '9HZC19EM0', '9HZC19LHG', '9HZC19LFG', '9HZC19LEG', '9HZC19EMI',
           'QIZC10EJ0', 'QIZC10EL0', 'QIZC10EM0', 'QIZC10LHG', 'QIZC10EMI',
           'QIZC11LHG', 'QIZC11LGG', 'QIZC11LFG', 'QIZC11LEG', 'QIZC11EMI',
           'QIZC12EI0', 'QIZC12EL0', 'QIZC12LHG', 'QIZC12LGG', 'QIZC12LFG',
           'QIZC12LEG', 'QIZC12EMI', 'QIZC13EI0', 'QIZC13EK0', 'QIZC13EL0',
           'QIZC13EJ0', 'QIZC13EM0', 'QIZC13LHG', 'QIZC13LFG', 'QIZC15EJ0',
           'QIZC13LEG', 'QIZC13EMI', 'QIZC14EMI', 'QIZC14EI0', 'QIZC15EI0',
           'QIZC14EK0', 'QIZC14EL0', 'QIZC14EM0', 'QIZC14LHG', 'QIZC14LFG',
           'QIZC14LGG', 'QIZC14LEG', 'QIZC15EK0', 'QIZC15EL0', 'QIZC15EMO',
           'QIZC15LGG', 'QIZC15LFG', 'QIZC16EL0', 'QIZC16EM0', 'QIZC15EMI',
           'QIZC16EI0', 'QIZC16LHG', 'QIZC16LGG', 'QIZC16LEG', 'QIZC16EMI',
           'QIZC17EI0', 'QIZC17EJ0', 'QIZC17EK0', 'QIZC17EL0', 'QIZC17EM0',
           'QIZC17LHG', 'QIZC17LGG', 'QIZC17LFG', 'QIZC18EJ0', 'QIZC17LEG',
           'QIZC17EMI', 'QIZC18EI0', 'QIZC18EK0', 'QIZC18LGG', 'QIZC18LFG',
           'QIZC18LEG', 'QIZC18EMI', 'QIZC19EI0', 'QIZC19EK0', 'QIZC19EL0',
           'QIZC19EM0', 'QIZC19LHG', 'QIZC19LGG', 'QIZC19LEG', 'QIZC19LFG',
           'QIZC19EJ0', 'QJZC10EI0', 'QIZC19EMI', 'QIZC15LEG', 'QLZC12EM0',
           'QIZC12EJ0', 'QIZC13LGG', 'QLZC23LGG', 'QLZC23LFG', 'QIZC16LFG',
           '9HHC13EL0', '9HHC13EM0', 'QIZC15LHG')

#PPID, age, gender, cohort assignment, race, ethnicity as columns
demo_df <- read_csv(file.path(".", "input_data", "demo", "STATE_Enrollment.csv"))
# cohort_sheets <- readxl::excel_sheets(file.path(".", "input_data", "cohorts", "STATE Draws - Deidentified_cohorts.xlsx"))
active <- read_excel(file.path(".", "input_data", "cohorts", "STATE Draws - Deidentified_cohorts.xlsx"), sheet = "Active")
active_ppids <- active$`Project Participant IDs`
past_cancer <- read_excel(file.path(".", "input_data", "cohorts", "STATE Draws - Deidentified_cohorts.xlsx"), sheet = "Past Cancer")
pc_ppids <- past_cancer$`Project Participant IDs`
nocancer <- read_excel(file.path(".", "input_data", "cohorts", "STATE Draws - Deidentified_cohorts.xlsx"), sheet = "No Cancer")
noc_ppids <- nocancer$`Project Participant IDs`

demo_df <- demo_df %>%
  mutate(white = recode(STATE_Race___White, "1" = "White", .default = ""),
         other = str_remove(STATE_Race_Other, "\\{null\\}"),
         afamer = recode(`STATE_Race___African American`, "1" = "African American", .default = ""),
         amerind = recode(`STATE_Race___American Indian or Alaska Native`, "1" = "American Indian or Alaska Native", .default = ""),
         asian = recode(STATE_Race___Asian, "1" = "Asian", .default = ""),
         dec = recode(`STATE_Race___Decline to Answer`, "1" = "Decline to Answer", .default = ""),
         namer = recode(`STATE_Race___Native American`, "1" = "Native American", .default = ""),
         ethnicity = recode(`STATE_Race___Hispanic or Latino`, "1" = "Hispanic or Latino", .default = "Not Hispanic or Latino")) %>%
  rowwise() %>%
  mutate(race = paste0(afamer, amerind, asian, namer, white, other, dec, collapse = ", ")) %>%
  ungroup() %>%
  mutate(race = case_when(race == "" ~ "Other",
                          TRUE ~ race),
         cohort = case_when(STATE_PPID %in% active_ppids ~ "Active",
                            STATE_PPID %in% pc_ppids ~ "Past Cancer",
                            STATE_PPID %in% noc_ppids ~ "No Cancer")) %>%
  dplyr::select(STATE_PPID, STATE_Age_Calculated, race, ethnicity, cohort) %>%
  filter(STATE_PPID %in% ppids)

write_csv(demo_df, file.path("input_data", "demographics.csv"))



