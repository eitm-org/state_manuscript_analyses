import os
import glob

from constants import workflow_results_dir, flat_results_dir

def main():
    hg002_fp = glob.glob(os.path.join(workflow_results_dir, '*', 'results', 'variant_calling/clairs', 'HG002*', 'fp/0000.vcf'))
    hg002_tp = glob.glob(os.path.join(workflow_results_dir, '*', 'results', 'variant_calling/clairs', 'HG002*', 'fp/0002.vcf'))

    if not os.path.exists(os.path.join(flat_results_dir, 'STATE_HG002_vcfs/fp')): 
        os.makedirs(os.path.join(flat_results_dir, 'STATE_HG002_vcfs/fp'))
    if not os.path.exists(os.path.join(flat_results_dir, 'STATE_HG002_vcfs/tp')): 
        os.makedirs(os.path.join(flat_results_dir, 'STATE_HG002_vcfs/tp'))
    for fp, tp in zip(hg002_fp, hg002_tp):
        fp = os.path.abspath(fp)
        tp = os.path.abspath(tp)
        base_name = fp.split('/')[-3].rsplit('_', 1)[0]
        fp_symlink_path = os.path.join(flat_results_dir, 'STATE_HG002_vcfs/fp', f'{base_name}.fp.vcf')
        tp_symlink_path = os.path.join(flat_results_dir, 'STATE_HG002_vcfs/tp', f'{base_name}.tp.vcf')
        if not os.path.islink(fp_symlink_path): 
            os.symlink(fp, fp_symlink_path)
        if not os.path.islink(tp_symlink_path): 
            os.symlink(tp, tp_symlink_path)

    hg002_mosdepth = glob.glob(os.path.join(workflow_results_dir, '*', 'results', 'reports/mosdepth', 'HG002*', '*.bed.gz'))
    if not os.path.exists(os.path.join(flat_results_dir, 'STATE_HG002_mosdepth')): 
        os.makedirs(os.path.join(flat_results_dir, 'STATE_HG002_mosdepth'))
    for mos in hg002_mosdepth:
        base_name = mos.split('/')[-1].rsplit('_', 1)[0]
        mos = os.path.abspath(mos)
        mos_symlink_path = os.path.join(flat_results_dir, 'STATE_HG002_mosdepth', f'{base_name}..md.regions.bed.gz')
        if not os.path.islink(mos_symlink_path): 
            os.symlink(mos, mos_symlink_path)

    eibs_mosdepth = glob.glob(os.path.join(workflow_results_dir, '*', 'results', 'reports/mosdepth', 'EIBS*', '*.bed.gz'))
    if not os.path.exists(os.path.join(flat_results_dir, 'STATE_mosdepth')): 
        os.makedirs(os.path.join(flat_results_dir, 'STATE_mosdepth'))
    for mos in eibs_mosdepth:
        file_name = mos.split('/')[-1]
        mos = os.path.abspath(mos)
        mos_symlink_path = os.path.join(flat_results_dir, 'STATE_mosdepth', file_name)
        if not os.path.islink(mos_symlink_path):
            os.symlink(mos, mos_symlink_path)

    eibs_vcf = glob.glob(os.path.join(workflow_results_dir, '*', 'results', 'variant_calling/clairs', 'EIBS*', '*.filtered3.vcf.gz'))
    if not os.path.exists(os.path.join(flat_results_dir, 'STATE_vcfs_f3/full_genome')): 
        os.makedirs(os.path.join(flat_results_dir, 'STATE_vcfs_f3/full_genome'))
    for vcf in eibs_vcf:
        file_name = vcf.split('/')[-1]
        vcf = os.path.abspath(vcf)
        vcf_symlink_path = os.path.join(flat_results_dir, 'STATE_vcfs_f3/full_genome', file_name)
        if not os.path.islink(vcf_symlink_path): 
            os.symlink(vcf, vcf_symlink_path)

    # eibs_vcf = glob.glob(os.path.join(workflow_results_dir, '*', 'results', 'variant_calling/clairs', 'EIBS*', '*.filtered1.vcf.gz'))
    # if not os.path.exists(os.path.join(flat_results_dir, 'STATE_vcfs_f1/full_genome')): 
    #     os.makedirs(os.path.join(flat_results_dir, 'STATE_vcfs_f1/full_genome'))
    # for vcf in eibs_vcf:
    #     file_name = vcf.split('/')[-1]
    #     vcf = os.path.abspath(vcf)
    #     vcf_symlink_path = os.path.join(flat_results_dir, 'STATE_vcfs_f1/full_genome', file_name)
    #     if not os.path.islink(vcf_symlink_path): 
    #         os.symlink(vcf, vcf_symlink_path)

if __name__ == "__main__":
    main()
