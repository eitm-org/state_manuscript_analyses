import glob
import os
import subprocess

from constants import refs_dir, flat_results_dir

def main():
    # NOTE: in your {refs_dir} download funcotator_dataSources.v1.7.20200521s from https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator/funcotator_dataSources.v1.7.20200521g?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))
    FUNC_CMD = """gatk Funcotator --variant {input_vcf} --reference {refs_dir}/GRCh38/GRCh38.primary_assembly.genome_X.fa \
        --ref-version hg38 --data-sources-path {refs_dir}/funcotator_dataSources.v1.7.20200521s --output {out_vcf} -output-file-format VCF
    """
    for patient_file in glob.glob(os.path.join(flat_results_dir, 'STATE_vcfs_f3_region_filtered/full_genome/*.vcf')):
        vcf = patient_file.split('/')[-1]
        vcf_basename, filternum = vcf.split('.')[0], vcf.split('.')[2]
        funcotated_dir = os.path.join(flat_results_dir, f'STATE_vcfs_f{filternum[-1]}_region_filtered_funcotated/full_genome/')
        if not os.path.exists(funcotated_dir):
            os.makedirs(funcotated_dir)
        output_path = os.path.join(funcotated_dir, vcf)
        sample = vcf_basename.split('_vs_')[0]
        if vcf in os.listdir(funcotated_dir):
            print(sample, 'skipped')
        else:
            print('Funcotating:', sample, vcf_basename, filternum)
            func_command = FUNC_CMD.format(input_vcf = patient_file, refs_dir = refs_dir, out_vcf = output_path)
            print(func_command)
            try:
                result = subprocess.check_output(func_command, shell=True, stderr=subprocess.STDOUT)
                print(result.decode('utf-8'))
            except subprocess.CalledProcessError as e:
                # error ino
                print(f"Error: {e.output.decode('utf-8')}")

if __name__ == "__main__":
    main()