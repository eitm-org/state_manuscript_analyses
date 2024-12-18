import glob
import os
import sys

import pandas as pd
import tensorsignatures as ts


from STATE_analyses.post_processing.utils import eid_to_patient_cohort, get_draw_mapping, parse_date
from aggregate_funcotator_utils import aggregate_funcotator_Mbps, aggregate_funcotator_chr, aggregate_funcotator_global, get_funcotator_data_for_vcf
from aggregate_vep_utils import aggregate_vep_Mbps, aggregate_vep_chrom, aggregate_vep_global, get_vep_df



def aggregate_vep(vep_global_path, vep_chrom_path, vep_Mbps_path, VEP_BASE_DIR):
    vep_dir = '/data/scratch/xchen/STATE_vcfs_f3_region_filtered/full_genome_vep/'
    vep_df = get_vep_df(vep_dir)
    vep_df_global = aggregate_vep_global(vep_df)
    vep_df_global.to_csv(vep_global_path)
    print(vep_df_global.shape)
    vep_df_chrom = aggregate_vep_chrom(vep_df)
    vep_df_chrom.to_csv(vep_chrom_path)
    vep_df_Mbps = aggregate_vep_Mbps(vep_df)
    vep_df_Mbps.to_csv(vep_Mbps_path)


def main():
    latest_run = str(sys.argv[1])
    filternums = ['f1', 'f3']
    for filternum in filternums:
        merged_global_patient_path = f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_snv_patient_global_{filternum}.csv'
        merged_chrom_patient_path = f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_snv_patient_chrom_{filternum}.csv'
        merged_Mbps_patient_path = f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_snv_patient_Mbps_{filternum}.csv'
        VCF_BASE_DIR = f'/data/scratch/xchen/STATE_vcfs_{filternum}_region_filtered_funcotated'
        VEP_BASE_DIR = f'/data/scratch/xchen/STATE_vcfs_{filternum}_region_filtered/full_genome_vep'

        funcotator_global_path = f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_funcotator_patient_global_{filternum}.csv'
        funcotator_chrom_path = f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_funcotator_patient_chrom_{filternum}.csv'
        funcotator_Mbps_path = f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_funcotator_patient_Mbps_{filternum}.csv'

        mutsig_global_path_prefix = f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/refit_state002_{latest_run}_{filternum}'
        mutsig_chrom_path_prefix = f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/refit_state002_{latest_run}_{filternum}_chrom'
        
        cosmic_global_path = f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/mut_sig_assignment_solution_activities_{filternum}.csv'
        
        kataegis_global_path = f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_kataegis_global_{filternum}.csv'
        kataegis_chrom_path = f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_kataegis_chrom_{filternum}.csv'
        kataegis_Mbps_path = f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_kataegis_Mbps_{filternum}.csv'

        vep_global_path = f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_vep_global_{filternum}.csv'
        vep_chrom_path = f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_vep_chrom_{filternum}.csv'
        vep_Mbps_path = f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_vep_Mbps_{filternum}.csv'

        # aggregate_funcotator(funcotator_global_path, funcotator_chrom_path, funcotator_Mbps_path, VCF_BASE_DIR)
        # aggregate_vep(vep_global_path, vep_chrom_path, vep_Mbps_path, VEP_BASE_DIR)
        aggregate_global(mutsig_global_path_prefix, funcotator_global_path, kataegis_global_path, vep_global_path, cosmic_global_path, merged_global_patient_path, VCF_BASE_DIR)
        aggregate_chrom(mutsig_chrom_path_prefix, funcotator_chrom_path, kataegis_chrom_path, vep_chrom_path, merged_chrom_patient_path, VCF_BASE_DIR)
        aggregate_Mbps(funcotator_Mbps_path, kataegis_Mbps_path, vep_Mbps_path, merged_Mbps_patient_path)



def aggregate_funcotator(funcotator_global_path, funcotator_chrom_path, funcotator_Mbps_path, VCF_BASE_DIR):
    vcf_full_genome_dir = os.path.join(VCF_BASE_DIR, 'full_genome')
    funcotator_global = aggregate_funcotator_global(vcf_full_genome_dir, debug=False)
    funcotator_global.to_csv(funcotator_global_path, index=False)
    funcotator_chr = aggregate_funcotator_chr(vcf_full_genome_dir, debug=False)
    funcotator_chr.to_csv(funcotator_chrom_path, index=False)
    funcotator_Mbps = aggregate_funcotator_Mbps(vcf_full_genome_dir, debug=False)
    funcotator_Mbps.to_csv(funcotator_Mbps_path, index=False)


def mutsig_global_to_df(mutsig_path_prefix, VCF_BASE_DIR) -> pd.DataFrame:
    input_vcf_paths2 = glob.glob(os.path.join(VCF_BASE_DIR, 'full_genome', '*.ann.filtered[1,3].vcf'))
    input_vcf_paths2.sort()
    patient_ids = [vcf_path.split('/')[-1].split('_vs_')[0] for vcf_path in input_vcf_paths2]
    i = 0
    exposure_concat_chunks = []
    for chunk in ['chunk1', 'chunk2', 'chunk3']:
        mutsig_path = f'{mutsig_path_prefix}_{chunk}.pkl'
        init = ts.load_dump(mutsig_path)
        rows = [f'ts{n+1:02}' for n in range(init.rank)]
        ids = patient_ids[i: i + init.sample_indices.shape[0]] #+ fp_ids + fp_pileup_ids
        exposure_df = pd.DataFrame(init.E.reshape(init.rank, init.sample_indices.shape[0]), columns = ids, index = rows).T
        exposure_df_norm = exposure_df.divide(exposure_df.sum(axis=1), axis=0)
        exposure_df = exposure_df.add_prefix('mutsig_count_')
        exposure_df_norm = exposure_df_norm.add_prefix('mutsig_ratio_')
        exposure_concat = pd.concat([exposure_df, exposure_df_norm], axis=1)
        exposure_concat['EIBS'] = exposure_concat.index
        exposure_concat_chunks.append(exposure_concat)
        i += init.sample_indices.shape[0]
    return pd.concat(exposure_concat_chunks, axis=0)

def mutsig_chrom_to_df(mutsig_path_prefix, VCF_BASE_DIR) -> pd.DataFrame:
    input_vcf_paths2 = glob.glob(os.path.join(VCF_BASE_DIR, 'chr*', '*.ann.filtered[1,3].chrom.vcf'))
    input_vcf_paths2.sort()
    patient_ids = [vcf_path.split('/')[-1].split('_vs_')[0] for vcf_path in input_vcf_paths2]
    chrs = [vcf_path.split('/')[-2] for vcf_path in input_vcf_paths2]
    i = 0
    exposure_concat_chunks = []
    for chunk in ['chunk1', 'chunk2', 'chunk3']:
        mutsig_path = f'{mutsig_path_prefix}_{chunk}.pkl'
        init = ts.load_dump(mutsig_path)
        rows = [f'ts{n+1:02}' for n in range(init.rank)]
        ids = patient_ids[i: i + init.sample_indices.shape[0]]
        chrs_chunk = chrs[i: i + init.sample_indices.shape[0]]
        exposure_df = pd.DataFrame(init.E.reshape(init.rank, init.sample_indices.shape[0]), columns = ids, index = rows).T
        exposure_df_norm = exposure_df.divide(exposure_df.sum(axis=1), axis=0)
        exposure_df = exposure_df.add_prefix('mutsig_count_')
        exposure_df_norm = exposure_df_norm.add_prefix('mutsig_ratio_')
        exposure_concat = pd.concat([exposure_df, exposure_df_norm], axis=1)
        exposure_concat['chr'] = chrs_chunk
        exposure_concat['EIBS'] = exposure_concat.index
        exposure_concat_chunks.append(exposure_concat)
        i += init.sample_indices.shape[0]
    return pd.concat(exposure_concat_chunks, axis=0)

def aggregate_global(
        mutsig_global_path,
        funcotator_global_path,
        kataegis_global_path,
        vep_global_path,
        cosmic_global_path,
        merged_global_patient_path,
        VCF_BASE_DIR
    ):
    exposure_global = mutsig_global_to_df(mutsig_global_path, VCF_BASE_DIR)
    funcotator_global = pd.read_csv(funcotator_global_path)
    kata_global = pd.read_csv(kataegis_global_path, index_col=0)
    kata_global = parse_date(kata_global)
    vep_global = pd.read_csv(vep_global_path)
    cosmic_global = pd.read_csv(cosmic_global_path)
    merged_global = (exposure_global
        .merge(kata_global, left_on=['EIBS'], right_on=['EIBS'], how='outer')
        .merge(funcotator_global, left_on=['EIBS'], right_on=['EIBS'], how='outer')
        .merge(vep_global, left_on=['EIBS'], right_on=['EIBS'], how='outer')
        .merge(cosmic_global, left_on=['EIBS'], right_on=['EIBS'], how='outer')
    )
    merged_global.index = merged_global.EIBS
    merged_global_patient = eid_to_patient_cohort(merged_global)
    merged_global_patient.reset_index(inplace=True)
    merged_global_patient['month'] = [6*dt if c == 'CONTROL' else 3*dt for dt, c in zip(merged_global_patient.draw_times, merged_global_patient.cohort)]
    merged_global_patient['month'] = merged_global_patient['month'].astype(str)
    merged_global_patient.to_csv(merged_global_patient_path, index=False)


def aggregate_chrom(
        mutsig_chrom_path,
        funcotator_chrom_path,
        kataegis_chrom_path,
        vep_chrom_path,
        merged_chrom_patient_path,
        VCF_BASE_DIR
):
    exposure_chrom = mutsig_chrom_to_df(mutsig_chrom_path, VCF_BASE_DIR)
    funcotator_chr = pd.read_csv(funcotator_chrom_path)
    kata_chrom = pd.read_csv(kataegis_chrom_path, index_col=0)
    kata_chrom = parse_date(kata_chrom)
    vep_chrom = pd.read_csv(vep_chrom_path)
    merged_chrom = (exposure_chrom
                    .merge(kata_chrom, left_on=['EIBS', 'chr'], right_on=['EIBS', 'chr'], how='outer')
                    .merge(funcotator_chr, left_on=['EIBS', 'chr'], right_on=['EIBS', 'chr'], how='outer')
                    .merge(vep_chrom, left_on=['EIBS', 'chr'], right_on=['EIBS', 'chr'], how='outer')
                )
    merged_chrom.index = merged_chrom.EIBS
    merged_chrom_patient = eid_to_patient_cohort(merged_chrom)
    merged_chrom_patient.reset_index(inplace=True)
    merged_chrom_patient['month'] = [6*dt if c == 'CONTROL' else 3*dt for dt, c in zip(merged_chrom_patient.draw_times, merged_chrom_patient.cohort)]
    merged_chrom_patient['month'] = merged_chrom_patient['month'].astype(str)
    merged_chrom_patient.to_csv(merged_chrom_patient_path, index=False)

def aggregate_Mbps(
        funcotator_Mbps_path,
        kataegis_Mbps_path,
        vep_Mbps_path,
        merged_Mbps_patient_path,
    ):
    funcotator_Mbps = pd.read_csv(funcotator_Mbps_path)
    kata_Mbps = pd.read_csv(kataegis_Mbps_path, index_col=0)
    kata_Mbps = parse_date(kata_Mbps)
    kata_Mbps['end'] = kata_Mbps.end.astype(int)
    vep_Mbps = pd.read_csv(vep_Mbps_path)
    merged_Mbps = (funcotator_Mbps
                    .merge(kata_Mbps, left_on=['EIBS', 'chr', 'start', 'end'], right_on=['EIBS', 'chr', 'start', 'end'], how='outer')
                    .merge(vep_Mbps, left_on=['EIBS', 'chr', 'start', 'end'], right_on=['EIBS', 'chr', 'start', 'end'], how='outer')
    )
    merged_Mbps.index = merged_Mbps['EIBS']
    merged_Mbps_patient = eid_to_patient_cohort(merged_Mbps)
    merged_Mbps_patient.reset_index(inplace=True)
    merged_Mbps_patient['month'] = [6*dt if c == 'CONTROL' else 3*dt for dt, c in zip(merged_Mbps_patient.draw_times, merged_Mbps_patient.cohort)]
    merged_Mbps_patient['month'] = merged_Mbps_patient['month'].astype(str)
    merged_Mbps_patient.to_csv(merged_Mbps_patient_path, index=False)


if __name__ == "__main__":
    main()
