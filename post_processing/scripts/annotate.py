import glob 
import os
import subprocess

from scripts.constants import workflow_results_dir, refs_dir

def main():
    FILTER_CMD = """bcftools view -f PASS {} -o {}\n"""
    FILL_STATS_CMD = """bcftools +fill-tags {} -Ob -o tmp.vcf.gz -- -t all,"INFO/GQ:1=int(sum(FORMAT/GQ))"\n"""
    MV_CMD = """mv tmp.vcf.gz {}\n"""
    INX_CMD = """bcftools index -t {} -f\n"""
    ANN_CMD = """
    gatk VariantAnnotator \
    -R /data/xchen/refs/GRCh38/GRCh38.primary_assembly.genome_X.fa \
    -V {input_path} \
    -O {output_path} \
    --resource:gnomad {refs_dir}/af-only-gnomad-remapped.vcf.gz \
    --resource:germline {germine_path}\
    -E gnomad.AF -E germline.F -E germline.P -E germline.GQ\
    --dbsnp {refs_dir}/broad_resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz
    """
    for patient_file in glob.glob(os.path.join(workflow_results_dir, 'STATE*/results/variant_calling/clairs/EIBS*/*_vs_*.clairs.vcf.gz')):
        vcf_dir, _ = patient_file.rsplit('/', 1)
        run, sample, vcf = patient_file.split('/')[-6], patient_file.split('/')[-2], patient_file.split('/')[-1]
        vcf_basename = vcf.split('.')[0]
        output_stats_dir = os.path.join(workflow_results_dir, run, 'results', 'variant_calling', 'clairs', sample)
        output_vcf_path  = os.path.join(output_stats_dir, f'{vcf_basename}.ann.vcf.gz')
        germline_path = os.path.join(output_stats_dir, f'{sample}_germline.clairs.vcf.gz')
        germline_pass_path = os.path.join(output_stats_dir, f'{sample}_germline.pass.clairs.vcf.gz')
        if  f'{vcf_basename}.ann.vcf.gz' in os.listdir(vcf_dir):
            print(f'{run} skipped')
            continue
        print(output_stats_dir, output_vcf_path, germline_path)
        filter_command = FILTER_CMD.format(germline_path, germline_pass_path)
        fill_stats_command = FILL_STATS_CMD.format(germline_pass_path)
        mv_command = MV_CMD.format(germline_pass_path)
        inx_command = INX_CMD.format(germline_pass_path)
        ann_command = ANN_CMD.format(
            input_path=patient_file, 
            output_path=output_vcf_path, 
            germine_path=germline_pass_path,
            refs_dir=refs_dir,
        )
        try:
            result = subprocess.check_output(filter_command, shell=True, stderr=subprocess.STDOUT)
            print(result.decode('utf-8'))
        except subprocess.CalledProcessError as e:
            # error ino
            print(f"Error: {e.output.decode('utf-8')}")
        try:
            result = subprocess.check_output(fill_stats_command, shell=True, stderr=subprocess.STDOUT)
            print(result.decode('utf-8'))
        except subprocess.CalledProcessError as e:
            # error ino
            print(f"Error: {e.output.decode('utf-8')}")
        try:
            result = subprocess.check_output(mv_command, shell=True, stderr=subprocess.STDOUT)
            print(result.decode('utf-8'))
        except subprocess.CalledProcessError as e:
            # error ino
            print(f"Error: {e.output.decode('utf-8')}")
        try:
            result = subprocess.check_output(inx_command, shell=True, stderr=subprocess.STDOUT)
            print(result.decode('utf-8'))
        except subprocess.CalledProcessError as e:
            # error ino
            print(f"Error: {e.output.decode('utf-8')}")
        try:
            result = subprocess.check_output(ann_command, shell=True, stderr=subprocess.STDOUT)
            print(result.decode('utf-8'))
        except subprocess.CalledProcessError as e:
            # error ino
            print(f"Error: {e.output.decode('utf-8')}")

if __name__ == "__main__":
    main()
