import glob 
import os
import subprocess

from constants import workflow_results_dir

def main():
    FILTER_CMD = """bcftools view -f PASS {} -o {}\n"""
    FILL_STATS_CMD = """bcftools +fill-tags {} -Ob -o tmp.vcf.gz -- -t all,"INFO/GQ:1=int(sum(FORMAT/GQ))"\n"""
    MV_CMD = """mv tmp.vcf.gz {}\n"""
    INX_CMD = """bcftools index -t {} -f\n"""
    ANN_CMD = """
    gatk VariantAnnotator \
    -R /data/xchen/refs/GRCh38/GRCh38.primary_assembly.genome_X.fa \
    -V {} \
    -O {} \
    --resource:gnomad /data/xchen/refs/af-only-gnomad-remapped.vcf.gz \
    --resource:germline {}\
    -E gnomad.AF -E germline.F -E germline.P -E germline.GQ\
    --dbsnp /data/scratch/xchen/STATE_analyses_data/broad_resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz
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
        ann_command = ANN_CMD.format(patient_file, output_vcf_path, germline_pass_path)
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


# VCFTOOLS_CMD = """
# vcftools --gzvcf {input} --out {output_prefix} {arg}
# """

# with open(os.path.join('backfill_vcftools_stats.sh'), 'w') as bash_file:
#     for hg002_file in glob.glob('/data/xchen/STATE_analyses_data/STATE*/results/variant_calling/clairs/HG002*/*germline.clairs.vcf.gz'):
#         print(hg002_file)
#         run, sample = hg002_file.split('/')[4], hg002_file.split('/')[-2]
#         output_stats_dir = os.path.join(base_path, run, 'results', 'reports', 'vcftools', 'clairs', sample)
#         make_result_dir(output_stats_dir)
#         bash_file.write(VCFTOOLS_CMD.format(input=hg002_file, output_prefix=os.path.join(output_stats_dir, f'{sample}.clairs'), arg="--TsTv-by-count"))
#         bash_file.write(VCFTOOLS_CMD.format(input=hg002_file, output_prefix=os.path.join(output_stats_dir, f'{sample}.clairs'), arg="--TsTv-by-qual"))
#         bash_file.write(VCFTOOLS_CMD.format(input=hg002_file, output_prefix=os.path.join(output_stats_dir, f'{sample}.clairs'), arg="--FILTER-summary"))
#         bash_file.write(VCFTOOLS_CMD.format(input=hg002_file, output_prefix=os.path.join(output_stats_dir, f'{sample}.clairs'), arg="--site-depth "))



# MULTIQC_CMD = "multiqc {input_dir} --config ~/projects/STATE_analyses/STATE_analyses/mulitqc_config.yml  --outdir {outdir} -f\n"
# with open(os.path.join('backfill_multiqc.sh'), 'w') as bash_file:
#     for result_dir in glob.glob('/data/xchen/STATE_analyses_data/STATE*/results'):
#         print(result_dir)
#         outdir = os.path.join(result_dir, 'multiqc')
#         bash_file.write(MULTIQC_CMD.format(input_dir=result_dir, outdir=outdir))

