
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if(!require(VariantAnnotation)){
    BiocManager::install("VariantAnnotation")
    library(VariantAnnotation)
}

if(!require(karyoploteR)){
     BiocManager::install("karyoploteR")
    library(karyoploteR)
}

if(!require(data.table)){
    install.packages("data.table")
    library(data.table)
}

if(!require(maftools)){
     BiocManager::install("maftools")
    library(maftools)
}

if(!require(Gmisc)){
    install.packages("Gmisc")
    library(Gmisc)
}

if(!require(reticulate)){
    install.packages("reticulate")
    library(reticulate)
}

source_python('constants.py')

eibs_vcf_path = file.path(flat_results_dir, '/STATE_vcfs_f1/full_genome')
region_filtered_vcf_path =  file.path(flat_results_dir, 'STATE_vcfs_f1_region_filtered/full_genome')

split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))

dir.create(file.path(flat_results_dir, 'STATE_vcfs_f1_region_filtered'), showWarnings = FALSE)
dir.create(file.path(flat_results_dir, 'STATE_vcfs_f1_region_filtered', 'full_genome'), showWarnings = FALSE)

blacklist_padding = 0
blacklist_path = file.path(bedfiles_dir, "problematic_regions/ENCODE_blacklist2.bed")
blacklist_data = fread(blacklist_path)
blacklist_data = blacklist_data[,1:4]
colnames(blacklist_data) =  c('chr', 'start', 'end', 'type')
blacklist_data$start = blacklist_data$start - blacklist_padding
blacklist_data$end = blacklist_data$end + blacklist_padding
blacklist = toGRanges(blacklist_data)



centromere_path = file.path(bedfiles_dir, "problematic_regions/centromeres.bed")
centromere_data = fread(centromere_path)
centromere_data = centromere_data[,1:4]
colnames(centromere_data) =  c('chr', 'start', 'end', 'type')
centromeres = toGRanges(centromere_data)

grc_path = file.path(bedfiles_dir, "problematic_regions/GRC_exclusions.bed")
grc_data = fread(grc_path)
grc_data = grc_data[,1:4]
colnames(grc_data) =  c('chr', 'start', 'end', 'type')
grc = toGRanges(grc_data)


unusual_path = file.path(bedfiles_dir, "problematic_regions/unusual_regions.bed")
unusual_data = fread(unusual_path)
unusual_data = unusual_data[,1:4]
colnames(unusual_data) =  c('chr', 'start', 'end', 'type')
unusual_grc_data = rbind(unusual_data, grc_data)
unusual_grc = toGRanges(unusual_grc_data)


repeats = file.path(bedfiles_dir, "problematic_regions/repeatmasker.bed")
repeats_data = fread(repeats)
repeats_data = repeats_data[,c(6,7,8,12)]
colnames(repeats_data) =  c('chr', 'start', 'end', 'type')
repeats_data = repeats_data[repeats_data$type == c('Satellite'),]
repeats = toGRanges(repeats_data)



adaptive_path = file.path(bedfiles_dir, 'adaptive/pathogenicGRCh38_20220906.bed')
adaptive_data = fread(adaptive_path)
adaptive_data = adaptive_data[,1:4]
colnames(adaptive_data) =  c('chr', 'start', 'end', 'type')
adaptive = toGRanges(adaptive_data)

subtract_data = rbind(unusual_grc_data, centromere_data, blacklist_data, adaptive_data)
subtract_ranges = toGRanges(subtract_data)


noisy_cov_basepath = file.path(flat_results_dir, 'noisy_coverage_regions')
eibs_vcf_files = list.files(eibs_vcf_path, pattern='*.vcf.gz$', full.names = T)
vcf_region_filtered_counts = c()
vcf_region_noise_filtered_counts = c()
eibs = c()
for(vcf_file in eibs_vcf_files){
  vcf_name = split_path(vcf_file)[1]
  vcf_name = strsplit(vcf_name, '.gz')[[1]][1]
  print(vcf_name)
  eibs_id = strsplit(vcf_name, '_vs_')[[1]][1]
  print(eibs_id)
  output_vcf_path  = paste0(region_filtered_vcf_path, '/', vcf_name)
  
  vcf <- readVcf(vcf_file)
  region_to_subtract = findOverlaps(vcf@rowRanges, subtract_ranges)
  vcf_region_filterd = vcf[-region_to_subtract@from]
  
  noisy_path = list.files(noisy_cov_basepath, pattern = eibs_id, full.names = T)[1]
  print(noisy_path)
  noisy_data = fread(noisy_path)
  noisy_data = noisy_data[,1:4]
  colnames(noisy_data) =  c('chr', 'start', 'end', 'type')
  noisy = toGRanges(noisy_data)
  noisy_to_subtract = findOverlaps(vcf_region_filterd@rowRanges, noisy)
  vcf_region_noise_filterd = vcf_region_filterd[-noisy_to_subtract@from]
  
  writeVcf(vcf_region_noise_filterd, output_vcf_path)
  
  print(paste('Before region+noise filter:', length(vcf)))
  print(paste('After region filter', length(vcf_region_filterd)))
  print(paste('After region+noise filter', length(vcf_region_noise_filterd)))
  vcf_region_filtered_counts = append(vcf_region_filtered_counts, length(vcf_region_filterd))
  vcf_region_noise_filtered_counts = append(vcf_region_noise_filtered_counts, length(vcf_region_noise_filterd))
  eibs = append(eibs, eibs_id)
}

region_filtered_counts = data.frame(eibs, vcf_region_filtered_counts, vcf_region_noise_filtered_counts)
write.csv(region_filtered_counts, 'region_noise_filtered_counts_f1.csv')
