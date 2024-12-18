
import glob
import os

from cyvcf2 import VCF, Writer

from STATE_analyses.post_processing.utils import ALL_CHRS


def split_vcf_by_chrom(input_vcf_path, output_vcf_paths):
    vcf = VCF(input_vcf_path)
    writers = [Writer(output_vcf_path, vcf) for output_vcf_path in output_vcf_paths]
    chrom_to_writers_map = dict(zip(ALL_CHRS, writers))
    for rec in vcf:
        writer = chrom_to_writers_map[rec.CHROM]
        writer.write_record(rec)
    for writer in chrom_to_writers_map.values(): writer.close()

def get_vcf_size(vcf):
    num_recs = 0
    for v in vcf: num_recs+=1
    return num_recs

def main():
    vcf_base_path = '/data/scratch/xchen/STATE_vcfs_f[1, 3]_region_filtered_funcotated'
    for input_vcf_path in glob.glob(os.path.join(vcf_base_path, 'full_genome', '*ann.filtered[1,3].vcf')):
        output_file_prefix = input_vcf_path.split('/')[-1][:-21]
        filternum = input_vcf_path.split('.')[2]
        input_base_path = input_vcf_path.split('full_genome')[0]
        output_file = output_file_prefix + f'.ann.{filternum}.chrom.vcf'
        output_file_dirs = [os.path.join(input_base_path, chr) for chr in ALL_CHRS]
        for dir in output_file_dirs:
            if not os.path.exists(dir): os.makedirs(dir)
        output_vcf_paths = [os.path.join(dir, output_file) for dir in output_file_dirs]
        split_vcf_by_chrom(input_vcf_path, output_vcf_paths)
    for input_vcf_path in glob.glob(os.path.join(vcf_base_path, 'full_genome', '*ann.filtered[1,3].vcf')):
        output_file_prefix = input_vcf_path.split('/')[-1][:-21]
        filternum = input_vcf_path.split('.')[2]
        input_base_path = input_vcf_path.split('full_genome')[0]
        output_file = output_file_prefix + f'.ann.{filternum}.chrom.vcf'
        output_file_dirs = [os.path.join(input_base_path, chr) for chr in ALL_CHRS]
        output_vcf_paths = [os.path.join(dir, output_file) for dir in output_file_dirs]
        for chrom_path in output_vcf_paths:
            vcf = VCF(chrom_path)
            if not get_vcf_size(vcf):
                print('Removing', chrom_path)
                os.remove(chrom_path)

if __name__ == "__main__":
    main()
