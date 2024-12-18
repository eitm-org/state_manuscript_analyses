import os
import glob

from cyvcf2 import VCF, Writer
import pandas as pd
from collections import defaultdict

from constants import workflow_results_dir

def apply_somatic_filter0(input_vcf_path, output_vcf_path, output_failed_vcf_path):
    mapping_filter_dict = {
        (False, False):  ('germline', 'quality'),
        (True, False): ('quality'),
        (False, True): ('germline'),
        (True, True): (),
    }
    filter_reasons_counts = defaultdict(int)
    vcf = VCF(input_vcf_path)
    w = Writer(output_vcf_path, vcf)
    w2 = Writer(output_failed_vcf_path, vcf)
    c = 0
    for rec in vcf:
        gnomad_af = rec.INFO.get('gnomad.AF') or -1
        germline_P = rec.INFO.get('germline.P') or False # False means it's not in germline
        germline_F = rec.INFO.get('germline.F') or False  # False means it's not in germline
        germline_filter = (gnomad_af < 0.0001) or not (germline_P or germline_F)
        quality_filter = rec.FILTER is None
        filter_reasons = mapping_filter_dict[(germline_filter, quality_filter)]
        filter_reasons_counts[filter_reasons]+=1
        if germline_filter and quality_filter:
            w.write_record(rec)
        else:
            w2.write_record(rec)
    w.close()
    w2.close()
    return filter_reasons_counts

def apply_somatic_filter1(input_vcf_path, output_vcf_path, output_failed_vcf_path, depth_chr_mapping):
    mapping_filter_dict = {
        (False, False, False):  ('germline', 'quality', 'depth'), # 000
        (False, False, True):  ('germline', 'quality'),           # 001

        (False, True, False): ('germline', 'depth'), # 010
        (False, True, True): ('germline'),           # 011

        (True, False, False): ('quality', 'depth'),  # 100
        (True, False, True): ('quality'),            # 101
        (True, True, False): ('depth'),              # 110                                            
        (True, True, True): (),                      # 111
    }
    filter_reasons_counts = defaultdict(int)
    vcf = VCF(input_vcf_path)
    w = Writer(output_vcf_path, vcf)
    w2 = Writer(output_failed_vcf_path, vcf)
    c = 0
    for rec in vcf:
        gnomad_af = rec.INFO.get('gnomad.AF') or -1
        germline_GQ = rec.INFO.get('germline.GQ') or -1
        germline_P = rec.INFO.get('germline.P') or False # False means it's not in germline
        germline_F = rec.INFO.get('germline.F') or False  # False means it's not in germline
        germline_F = (germline_GQ > 15) and germline_F
        germline_filter = (gnomad_af < 0.0001) and not (germline_P or germline_F)
        quality_filter = rec.FILTER is None
        depth_filter = rec.format('DP').sum() < (1.5 * depth_chr_mapping[rec.CHROM])
        filter_reasons = mapping_filter_dict[(germline_filter, quality_filter, depth_filter)]
        filter_reasons_counts[filter_reasons]+=1
        if germline_filter and quality_filter and depth_filter:
            w.write_record(rec)
        else:
            w2.write_record(rec)
    w.close()
    w2.close()
    return filter_reasons_counts

def apply_somatic_filter2(input_vcf_path, output_vcf_path, output_failed_vcf_path):
    mapping_filter_dict = {
        (False, False):  ('germline', 'quality'),
        (True, False): ('quality'),
        (False, True): ('germline'),
        (True, True): (),
    }
    filter_reasons_counts = defaultdict(int)
    vcf = VCF(input_vcf_path)
    w = Writer(output_vcf_path, vcf)
    w2 = Writer(output_failed_vcf_path, vcf)
    c = 0
    for rec in vcf:
        gnomad_af = rec.INFO.get('gnomad.AF') or -1
        germline_GQ = rec.INFO.get('germline.GQ') or -1
        germline_P = rec.INFO.get('germline.P') or False # False means it's not in germline
        germline_F = rec.INFO.get('germline.F') or False  # False means it's not in germline
        germline_F = (germline_GQ > 15) and germline_F
        germline_filter = (gnomad_af < 0) and not (germline_P or germline_F)
        quality_filter = rec.FILTER is None
        filter_reasons = mapping_filter_dict[(germline_filter, quality_filter)]
        filter_reasons_counts[filter_reasons]+=1
        if germline_filter and quality_filter:
            w.write_record(rec)
        else:
            w2.write_record(rec)
    w.close()
    w2.close()
    return filter_reasons_counts

def apply_somatic_filter3(input_vcf_path, output_vcf_path, output_failed_vcf_path, depth_chr_mapping):
    mapping_filter_dict = {
        (False, False, False):  ('germline', 'quality', 'depth'), # 000
        (False, False, True):  ('germline', 'quality'),           # 001

        (False, True, False): ('germline', 'depth'), # 010
        (False, True, True): ('germline'),           # 011

        (True, False, False): ('quality', 'depth'),  # 100
        (True, False, True): ('quality'),            # 101
        (True, True, False): ('depth'),              # 110                                            
        (True, True, True): (),                      # 111
    }
    filter_reasons_counts = defaultdict(int)
    vcf = VCF(input_vcf_path)
    w = Writer(output_vcf_path, vcf)
    w2 = Writer(output_failed_vcf_path, vcf)
    c = 0
    for rec in vcf:
        gnomad_af = rec.INFO.get('gnomad.AF') or -1
        germline_P = rec.INFO.get('germline.P') or False # False means it's not in germline
        germline_F = rec.INFO.get('germline.F') or False  # False means it's not in germline
        germline_filter = (gnomad_af < 0) and not (germline_P or germline_F)
        quality_filter = rec.FILTER is None

        depth_filter = rec.format('DP').sum() < (1.5 * depth_chr_mapping[rec.CHROM])
        filter_reasons = mapping_filter_dict[(germline_filter, quality_filter, depth_filter)]
        filter_reasons_counts[filter_reasons]+=1
        if germline_filter and quality_filter and depth_filter:
            w.write_record(rec)
        else:
            # rec.INFO['filter_reason'] = ((gnomad_af < 0) and not (germline_P or germline_F)
            w2.write_record(rec)
    w.close()
    w2.close()
    return filter_reasons_counts



def main():
    filter1_reasons_counts_all_samples = {}
    filter3_reasons_counts_all_samples = {}

    for input_vcf_path in glob.glob(os.path.join(workflow_results_dir, '*', 'results/variant_calling/clairs', '*', '*.ann.vcf.gz')):
        vcf_dir, vcf_filename = input_vcf_path.rsplit('/', 1)
        # if glob.glob(os.path.join(vcf_dir, '*.filtered1.vcf.gz')) and glob.glob(os.path.join(vcf_dir, '*.filtered3.vcf.gz')): continue
        vcf_basename = vcf_filename.split('.')[0]
        eibs = vcf_basename.split('_vs_')[0]
        run_id = vcf_dir.split('/')[-5]
        print(run_id, vcf_basename)
        mosdepth_path = os.path.join(workflow_results_dir, run_id, 'results/reports/mosdepth', eibs, f'{eibs}.md.mosdepth.summary.txt')
        if not os.path.exists(mosdepth_path): # in case some mosdepth files don't exist
            print(mosdepth_path)
            break
        depth = pd.read_csv(mosdepth_path, sep='\t')
        depth_chr = depth.loc[~depth.chrom.str.endswith('region') & depth.chrom.str.startswith('chr')]
        depth_chr_mapping = {chr:depth for chr, depth in zip(depth_chr.chrom, depth['mean'])}

        output_vcf_path = os.path.join(vcf_dir, f'{vcf_basename}.ann.filtered1.vcf.gz')
        output_failed_vcf_path = os.path.join(vcf_dir, f'{vcf_basename}.ann.failed_filter1.vcf.gz')
        # if (output_vcf_path in os.listdir(vcf_dir)) and (output_failed_vcf_path in os.listdir(vcf_dir)): continue
        filter_reasons_counts = apply_somatic_filter1(input_vcf_path, output_vcf_path, output_failed_vcf_path, depth_chr_mapping)
        filter1_reasons_counts_all_samples[vcf_basename] = filter_reasons_counts
        print('f1', filter_reasons_counts)

        output_vcf_path = os.path.join(vcf_dir, f'{vcf_basename}.ann.filtered3.vcf.gz')
        output_failed_vcf_path = os.path.join(vcf_dir, f'{vcf_basename}.ann.failed_filter3.vcf.gz')
        # if (output_vcf_path in os.listdir(vcf_dir)) and (output_failed_vcf_path in os.listdir(vcf_dir)): continue
        filter_reasons_counts = apply_somatic_filter3(input_vcf_path, output_vcf_path, output_failed_vcf_path, depth_chr_mapping)
        filter3_reasons_counts_all_samples[vcf_basename] = filter_reasons_counts
        print('f3', filter_reasons_counts)
        print(output_vcf_path)
    pd.DataFrame.from_dict(filter1_reasons_counts_all_samples).T.to_csv('filter1_reasons_counts_all_samples.csv', mode='a')
    pd.DataFrame.from_dict(filter3_reasons_counts_all_samples).T.to_csv('filter3_reasons_counts_all_samples.csv', mode='a')


if __name__ == "__main__":
    main()