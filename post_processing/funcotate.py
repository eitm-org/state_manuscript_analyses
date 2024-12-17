import glob
import os
import subprocess


def main():
    FUNC_CMD = """gatk Funcotator --variant {} --reference /data/xchen/refs/GRCh38/GRCh38.primary_assembly.genome_X.fa \
        --ref-version hg38 --data-sources-path /data/xchen/refs/funcotator_dataSources.v1.7.20200521s --output {} -output-file-format VCF
    """

    for patient_file in glob.glob('/data/scratch/xchen/STATE_vcfs_f[1,3]_region_filtered/full_genome/*.vcf'):
        vcf_dir, _ = patient_file.rsplit('/', 1)
        vcf = patient_file.split('/')[-1]
        vcf_basename, filternum = vcf.split('.')[0], vcf.split('.')[2]
        funcotated_dir = f'/data/scratch/xchen/STATE_vcfs_f{filternum[-1]}_region_filtered_funcotated/full_genome/'
        output_path = os.path.join(funcotated_dir, vcf)
        sample = vcf_basename.split('_vs_')[0]
        # print(os.listdir(funcotated_dir))
        if vcf in os.listdir(funcotated_dir):
            print(sample, 'skipped')
        else:
            print('Funcotating:', sample, vcf_basename, filternum)
            func_command = FUNC_CMD.format(patient_file, output_path)
            print(func_command)
            try:
                result = subprocess.check_output(func_command, shell=True, stderr=subprocess.STDOUT)
                print(result.decode('utf-8'))
            except subprocess.CalledProcessError as e:
                # error ino
                print(f"Error: {e.output.decode('utf-8')}")

if __name__ == "__main__":
    main()