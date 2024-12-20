import glob
import os
import subprocess

from scripts.constants import HG002_fully_resolved_path, workflow_results_dir

def main():
    ISEC_CMD = """bcftools isec -f PASS {input_file} {HG002_fully_resolved_path} -p {output_path}\n"""
    for input_file in glob.glob(os.path.join(workflow_results_dir, '*', 'results/variant_calling/clairs', 'HG002_*', '*_germline.clairs.vcf.gz')):
        vcf_dir, vcf_filename = input_file.rsplit('/', 1)
        if 'fp' in os.listdir(vcf_dir):
            continue
        vcf_basename = vcf_filename.split('.')[0]
        isec_cmd = ISEC_CMD.format(input_file = input_file, HG002_fully_resolved_path = HG002_fully_resolved_path, output_path = os.path.join(vcf_dir, 'fp'))
        try:
            result = subprocess.check_output(isec_cmd, shell=True, stderr=subprocess.STDOUT)
            print(result.decode('utf-8'))
        except subprocess.CalledProcessError as e:
            # error ino
            print(f"Error: {e.output.decode('utf-8')}")
        print(f'Output path: {os.path.join(vcf_dir, "fp")}')
if __name__ == "__main__":
    main()
