import os
import glob
import subprocess

from scripts.constants import flat_results_dir

def main():
    VEP_COMMAND = "vep -i {} -o {} --cache --force_overwrite --sift b --polyphen b -v --plugin AlphaMissense,file=/data/xchen/refs/vep/AlphaMissense_hg38.tsv.gz,cols=all"
    filternums = ['f3']
    for filternum in filternums:
        input_dir = os.path.join(flat_results_dir, f'STATE_vcfs_{filternum}_region_filtered/full_genome')
        output_dir = os.path.join(flat_results_dir, f'STATE_vcfs_{filternum}_region_filtered/full_genome_vep')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for input_file in glob.glob(os.path.join(input_dir, '*.vcf')):
            vcf_dir, vcf_filename = input_file.rsplit('/', 1)
            output_file = os.path.join(output_dir, vcf_filename)
            if os.path.exists(output_file):
                continue
            vep_cmd = VEP_COMMAND.format(input_file, output_file)
            print(vep_cmd)
            try:
                result = subprocess.check_output(vep_cmd, shell=True, stderr=subprocess.STDOUT)
                print(result.decode('utf-8'))
            except subprocess.CalledProcessError as e:
                # error ino
                print(f"Error: {e.output.decode('utf-8')}")

if __name__ == "__main__":
    main()