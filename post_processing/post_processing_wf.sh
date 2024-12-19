#!/usr/bin/env bash
echo "extract false/true positives from HG002"
python scripts/get_hg002_fp.py
echo "annotate patient somatic vcfs with germline, dbsnp, etc"
python scripts/annotate.py
echo "apply somatic filter"
python scripts/apply_somatic_filter.py
echo "symlink all relevant files to flat directory structure"
python scripts/collect_files.py
echo "apply noise model to save noisy coverage bedfiles for filtration"
python scripts/get_noisy_coverage_filter.py
echo "apply region-based and noisy coverage-based filtration"
Rscript scripts/region_noise_filter_f3.R
echo "add functional annotations"
python scripts/funcotate.py
echo "refit tensorignatures"
python scripts/refit_tensig.py
python scripts/get_vep.py
echo "assign cosmic signatures"
python scripts/assign_cosmic.py 
echo "aggregate SNV metrics"
python scripts/aggregate_snv_metrics.py
echo "batch correct"
Rscript scripts/batch_effect_correction_functions_f3.R
