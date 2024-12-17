#!/usr/bin/env bash

echo "rclone pipeline results from oci"
rclone copy oci:Results_Restricted/ --include **/mosdepth/** --include **/variant_calling/** --include **modkit** /data/scratch/xchen/STATE_analyses_data --verbose
echo "extract false/true positives from HG002"
python get_hg002_fp.py
echo "annotate patient somatic vcfs with germline, dbsnp, etc"
python annotate.py
echo "apply somatic filter"
python apply_somatic_filter.py
echo "symlink all relevant files to flat directory structure"
python collect_files.py
echo "apply noise model to save noisy coverage bedfiles for filtration"
python get_noisy_coverage_filter.py
echo "apply region-based and noisy coverage-based filtration"
Rscript region_noise_filter_f3.R
Rscript region_noise_filter_f1.R
echo "add functional annotations" # TODO: fix the rest 
python ~/projects/STATE_analyses/STATE_analyses/post_processing/funcotate.py
echo "split vcfs by chromosome"
python ~/projects/STATE_analyses/STATE_analyses/post_processing/split_vcfs_by_chrom.py
echo "calculate and aggregate kataegis metrics"
Rscript ~/projects/STATE_analyses/R/hg002_fp_kataegis.R
Rscript ~/projects/STATE_analyses/R/aggregate_katagis_global.R
Rscript ~/projects/STATE_analyses/R/aggregate_kataegis_chrom.R
Rscript ~/projects/STATE_analyses/R/aggregate_katagis_Mbps.R
echo "refit tensorignatures"
python ~/projects/STATE_analyses/STATE_analyses/post_processing/refit_tensig.py $1
python  ~/projects/STATE_analyses/STATE_analyses/post_processing/get_vep.py
echo "assign cosmic signatures"
python ~/projects/STATE_analyses/STATE_analyses/post_processing/assign_cosmic.py 
echo "aggregate SNV metrics"
python ~/projects/STATE_analyses/STATE_analyses/post_processing/aggregate_snv_metrics.py $1
echo "batch correct"
Rscript ~/projects/STATE_QC/RScripts/batch_effect_correction_functions_f3.R
Rscript ~/projects/STATE_QC/RScripts/batch_effect_correction_functions_f1.R
