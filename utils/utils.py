import os
from collections import defaultdict

from cyvcf2 import VCF, Writer
import glob
import pandas as pd
import numpy as np
import scipy

import tensorsignatures as ts

repeats = ['EIBS-001US', 'EIBS-001TB', 'EIBS-001TA', 'EIBS-001YM']

BASELINE = 0
DRAW2 = 1
DRAW3 = 2
DRAW4 = 3
DRAW5 = 4
DRAW6 = 5

COHORT_WELLNESS = 'WELLNESS'
COHORT_STAFF = 'STAFF'

COHORT_TMT = 'TREATMENT'
COHORT_CHEMO = 'CHEMO'
COHORT_MAPPING = {'Treatment': COHORT_TMT, 'Nonpatient': COHORT_STAFF, 'Wellness': COHORT_WELLNESS, 'Chemotherapy': COHORT_CHEMO}
ALL_CHRS = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']


COHORT_NC = 'no_cancer'
COHORT_PC = 'past_cancer'
COHORT_AC = 'active_cancer'
CLINICAL_COHORT_MAPPING = {'No Cancer': COHORT_NC, 'Past Cancer': COHORT_PC, 'Active': COHORT_AC}


aggregate_snv_Mbps_f3_path = '/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_snv_patient_Mbps_f3.csv'
aggregate_snv_chrom_f3_path = '/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_snv_patient_chrom_f3.csv'
aggregate_snv_global_f3_path = '/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_snv_patient_global_f3.csv'

aggregate_snv_Mbps_f1_path = '/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_snv_patient_Mbps_f1.csv'
aggregate_snv_chrom_f1_path = '/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_snv_patient_chrom_f1.csv'
aggregate_snv_global_f1_path = '/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/aggregate_snv_patient_global_f1.csv'
sbs_artifacts = ['SBS27', 'SBS43', 'SBS45', 'SBS46', 'SBS47', 'SBS48', 'SBS49', 'SBS50', 'SBS51', 'SBS52', 'SBS53', 'SBS54', 'SBS55', 'SBS56', 'SBS57', 'SBS58', 'SBS59', 'SBS60', 'SBS95']

def derive_final_metrics_uncorrrected(data, error_model=False):
    data = data.drop_duplicates()
    sbs_cols = [col for col in data.columns if col.startswith('SBS') and not col.endswith('adj3') and col.split('_')[0] not in sbs_artifacts]
    if error_model:
        error_model = data[[ s + '_count' for s in sbs_artifacts]].sum(axis=1) / data[sbs_cols].sum(axis=1)
        data['error_model'] = error_model
        count_cols = [col for col in data.columns if 'count' in col and 'SBS' not in col]
        data[count_cols] = data[count_cols].multiply(data.error_model, axis='index')
    mutsig_cols = [col for col in data.columns if col.startswith('mutsig') and not col.endswith('adj3')]
    
    sbs_norm = data[sbs_cols].divide(data[sbs_cols].sum(axis=1), axis=0)
    mutsig_norm = data[mutsig_cols].divide(data[mutsig_cols].sum(axis=1), axis=0)
    sbs_norm.columns = sbs_norm.columns.str.replace('count', 'ratio')
    mutsig_norm.columns = mutsig_norm.columns.str.replace('count', 'ratio')
    sbs_norm_cols = sbs_norm.columns.tolist()
    mutsig_norm_cols = mutsig_norm.columns.tolist()

    other_adj3_count_cols = [col for col in data.columns if 'count' in col and not col.endswith('adj3') and 'mutsig' not in col and not col.startswith('SBS')]
    other_adj3_count_metrcs = data[other_adj3_count_cols]
    other_adj3_ratio_metrics = other_adj3_count_metrcs.divide(data.snv_count, axis=0)
    other_adj3_ratio_metrics.columns = other_adj3_ratio_metrics.columns.str.replace('count', 'ratio')
    data = pd.concat([data, sbs_norm, mutsig_norm, other_adj3_ratio_metrics], axis=1)
    num = data._get_numeric_data()
    num[num < 0] = 0
    return data


def derive_final_metrics(data, error_model=True):
    data = data.drop_duplicates()
    sbs_cols = [col for col in data.columns if col.startswith('SBS') and col.endswith('adj3') and col.split('_')[0] not in sbs_artifacts]
    if error_model:
        error_model = data[[ s + '_count_adj3' for s in sbs_artifacts]].sum(axis=1) / data[sbs_cols].sum(axis=1)
        data['error_model'] = error_model
        count_cols = [col for col in data.columns if 'count' in col and 'SBS' not in col]
        data[count_cols] = data[count_cols].multiply(data.error_model, axis='index')
    mutsig_cols = [col for col in data.columns if col.startswith('mutsig') and col.endswith('adj3')]
    
    sbs_norm = data[sbs_cols].divide(data[sbs_cols].sum(axis=1), axis=0)
    mutsig_norm = data[mutsig_cols].divide(data[mutsig_cols].sum(axis=1), axis=0)
    sbs_norm.columns = sbs_norm.columns.str.replace('count', 'ratio')
    mutsig_norm.columns = mutsig_norm.columns.str.replace('count', 'ratio')
    sbs_norm_cols = sbs_norm.columns.tolist()
    mutsig_norm_cols = mutsig_norm.columns.tolist()

    other_adj3_count_cols = [col for col in data.columns if 'count' in col and col.endswith('adj3') and 'mutsig' not in col and not col.startswith('SBS')]
    other_adj3_count_metrcs = data[other_adj3_count_cols]
    other_adj3_ratio_metrics = other_adj3_count_metrcs.divide(data.snv_count_adj3, axis=0)
    other_adj3_ratio_metrics.columns = other_adj3_ratio_metrics.columns.str.replace('count', 'ratio')
    data = pd.concat([data, sbs_norm, mutsig_norm, other_adj3_ratio_metrics], axis=1)
    num = data._get_numeric_data()
    num[num < 0] = 0
    return data

def read_agg_snv(path):
    agg_snv = pd.read_csv(path)
    agg_snv['date'] = pd.to_datetime(agg_snv.date, format='%Y-%m-%d', errors='ignore')
    agg_snv['date'] = agg_snv['date'].dt.date
    # agg_snv.sort_values(['month', 'date'], inplace=True)
    if 'chr' in agg_snv.columns:
        agg_snv['chr_num'] = pd.to_numeric(agg_snv.chr.str[3:].replace('X', '23').replace('Y', '24'))
    agg_snv.index=agg_snv.EIBS
    for col in ['draw_id', 'draw_times', 'cohort']:
        if col in agg_snv.columns:
            agg_snv.drop(columns=col, inplace=True)
    return eid_to_patient_cohort(agg_snv).reset_index()

def join_clinical(agg_snv):
    # rcc_path = '/home/xchen@okta-oci.eitm.org/dropbox/Redcap_Cloud_STATE/rcc_user_records_20240124_0923/STATE_Enrollment/STATE_Enrollment.csv'
    # rcc = pd.read_csv(rcc_path)
    # agg_snv = agg_snv.merge(rcc, left_on='draw_id', right_on='participantId', how='left')
    agg_snv['age_high'] = agg_snv.draw_age >= 65
    htx_path = '/home/xchen@okta-oci.eitm.org/dropbox/Redcap_Cloud_STATE/archive/rcc_user_records_20240603_0922/STATE_History_and_Diagnosis/STATE_History.csv'
    htx = pd.read_csv(htx_path)[['STATE_Alc_Consumption', 'STATE_Tobacco_History', 'participantId']]
    br_htx_path = '/home/xchen@okta-oci.eitm.org/dropbox/Redcap_Cloud_STATE/rcc_user_records_20240124_0923/STATE_History_and_Diagnosis_Breast/STATE_History.csv'
    br_htx = pd.read_csv(br_htx_path)[['STATE_Alc_Consumption', 'STATE_Tobacco_History', 'participantId']]
    htx = pd.concat([htx, br_htx], axis=0)
    agg_snv = agg_snv.merge(htx, left_on='draw_id', right_on='participantId', how='left')
    # agg_snv.drop(columns='cohort_other')
    return agg_snv

def get_select_draw_ids():
    draw_mapping = get_draw_mapping()
    select_draw_ids = []
    for key, (did, dt, co) in draw_mapping.items():
        if dt >= 3 and did not in select_draw_ids:
            select_draw_ids.append(did)
    return select_draw_ids

def get_draw_mapping():
    mapping_path = '/home/xchen@okta-oci.eitm.org/dropbox/STATE_Draw_Event_Data/STATE Draws - Deidentified_cohorts.xlsx'
    draw_mapping = {}
    cohort_sheets = ['No Cancer', 'Past Cancer', 'Active']
    for sheet in cohort_sheets:
        draw_table = pd.read_excel(mapping_path, sheet_name=sheet)
        for i, row in draw_table.iterrows():
            
            patient_id = row['Project Participant IDs']
            age = row['Age']
            baseline_eid = row['Baseline']
            draw2_eid = row['Draw 2']
            draw3_eid = row['Draw 3']
            draw4_eid = row['Draw 4']
            draw5_eid = row['Draw 5']
            draw6_eid = row['Draw 6']
            draw1_mon = 0
            draw2_mon = row['Month since Baseline \n(Draw 2)']
            draw3_mon = row['Month since Baseline \n(Draw 3)']
            draw4_mon = row['Month since Baseline \n(Draw 4)']
            draw5_mon = row['Month since Baseline \n(Draw 5)']
            draw6_mon = row['Month since Baseline \n(Draw 6)']
            # cohort = COHORT_MAPPING[row['Cohorts']]
            cohort = CLINICAL_COHORT_MAPPING[sheet]
            draw_mapping[baseline_eid] = (patient_id, BASELINE, cohort, draw1_mon, age)
            draw_mapping[baseline_eid] = (patient_id, BASELINE, cohort, draw1_mon, age)
            if draw2_eid is not np.nan:
                draw_mapping[draw2_eid] = (patient_id, DRAW2, cohort, draw2_mon, age)
                draw_mapping[draw2_eid] = (patient_id, DRAW2, cohort, draw2_mon, age)
            if draw3_eid is not np.nan:
                draw_mapping[draw3_eid] = (patient_id, DRAW3, cohort, draw3_mon, age)
                draw_mapping[draw3_eid] = (patient_id, DRAW3, cohort, draw3_mon, age)
            if draw4_eid is not np.nan:
                draw_mapping[draw4_eid] = (patient_id, DRAW4, cohort, draw4_mon, age)
                draw_mapping[draw4_eid] = (patient_id, DRAW4, cohort, draw4_mon, age)
            if draw5_eid is not np.nan:
                draw_mapping[draw5_eid] = (patient_id, DRAW5, cohort, draw5_mon, age)
                draw_mapping[draw5_eid] = (patient_id, DRAW5, cohort, draw5_mon, age)
            if draw6_eid is not np.nan:
                draw_mapping[draw6_eid] = (patient_id, DRAW6, cohort, draw6_mon, age)
                draw_mapping[draw6_eid] = (patient_id, DRAW6, cohort, draw6_mon, age)
    return draw_mapping


def eid_to_patient_cohort(df):
    draw_mapping = get_draw_mapping()
    df_copy = df.copy()
    draw_id = [draw_mapping[idx[0]][0] for idx in df_copy.index.str.split('_', n=1)]
    draw_times = [draw_mapping[idx[0]][1] for idx in df_copy.index.str.split('_', n=1)]
    cohorts = [draw_mapping[idx[0]][2] for idx in df_copy.index.str.split('_', n=1)]
    draw_month = [draw_mapping[idx[0]][3] for idx in df_copy.index.str.split('_', n=1)]
    draw_age = [draw_mapping[idx[0]][4] for idx in df_copy.index.str.split('_', n=1)]
    df_copy.index = pd.MultiIndex.from_tuples(zip(draw_id, draw_times, cohorts, draw_month, draw_age), names=["draw_id", "draw_times", "cohort", "draw_month", "draw_age"])
    # df_copy.sort_index(inplace=True)
    return df_copy

def eid_to_repeats(df):
    df_copy = df.copy()
    repeats = ['EIBS-001US', 'EIBS-001TB', 'EIBS-001TA']
    df_copy.index = [idx if idx in repeats else 'other' for idx in df_copy.index.str.split('_', expand=True).get_level_values(0)]
    df_copy.sort_index(inplace=True)
    return df_copy


def parse_date(df):
    dates = (df\
    .index.str.replace('_[0-9]_', '_', regex=True)\
    .str.replace('_[0-9][0-9]_', '_', regex=True)\
    .str.split('_', expand=True).get_level_values(3))
    df['date'] = pd.to_datetime(dates, format='%Y%m%d', errors='ignore')
    df['date'] = df['date'].dt.date
    return df.sort_values('date')


def read_vcf_to_pandas(filepath, contig=None):
    vcf_dict = defaultdict(list)
    vcf = VCF(filepath)
    if contig is not None:
        vcf = vcf(contig)
    for v in vcf:
        update_vcf_dict(vcf_dict, v)
    return pd.DataFrame(vcf_dict)

def update_vcf_dict(vcf_dict, v):
    GT_mapping = {
        '\x02\x04': '0/1',
        '\x04\x04': '1/1',
        '\x02\x02': '0/0',
        '\x04\x02': '1/0',
    }
    vcf_dict['CHROM'].append(v.CHROM)
    vcf_dict['POS'].append(v.POS)
    vcf_dict['REF'].append(v.REF)
    vcf_dict['ALT'].append(v.ALT[0])
    vcf_dict['QUAL'].append(v.QUAL)
    vcf_dict['DBSNP_ID'].append(v.ID)
    vcf_dict['gnomad.AF'].append(v.INFO.get('gnomad.AF') or -1)
    vcf_dict['germline_P'].append(v.INFO.get('germline.P') or False)
    vcf_dict['germline_F'].append(v.INFO.get('germline.F') or False)
    vcf_dict['germline.GQ'].append(v.INFO.get('germline.GQ') or -1)
    vcf_dict['GT'].append(GT_mapping[v.format('GT')[0]])
    vcf_dict['DP'].append(v.format('DP').sum())
    vcf_dict['GQ'].append(v.format('GQ').sum())
    vcf_dict['AF'].append(v.format('AF').sum())
    vcf_dict['MAF'].append(v.INFO.get('MAF'))
    vcf_dict['ExcHet'].append(v.INFO.get('ExcHet'))
    vcf_dict['cosmic.COSMIC_ID'].append(v.INFO.get('cosmic.COSMIC_ID'))
    vcf_dict['cosmic.COSMIC_GENE'].append(v.INFO.get('cosmic.COSMIC_GENE'))
    vcf_dict['cosmic.COSMIC_STRAND'].append(v.INFO.get('cosmic.COSMIC_STRAND'))
    vcf_dict['AU'].append(v.format('AU').sum())
    vcf_dict['CU'].append(v.format('CU').sum())
    vcf_dict['GU'].append(v.format('GU').sum())
    vcf_dict['TU'].append(v.format('TU').sum())
    vcf_dict['NDP'].append(v.format('NDP').sum())
    vcf_dict['NAF'].append(v.format('NAF').sum())
    vcf_dict['FILTER'].append(v.FILTER)

def update_vcf_germline_dict(vcf_dict, v):
    vcf_dict['CHROM'].append(v.CHROM)
    vcf_dict['POS'].append(v.POS)
    vcf_dict['QUAL'].append(v.QUAL)
    vcf_dict['F'].append(v.INFO.get('F') or False)
    vcf_dict['P'].append(v.INFO.get('P') or False)
    vcf_dict['DP'].append(v.format('DP').sum())
    vcf_dict['GQ'].append(v.format('GQ').sum())
    vcf_dict['AF'].append(v.format('AF').sum())
    vcf_dict['FILTER'].append(v.FILTER)


def calculate_slopes(agg_snv_global, cols, meta_cols):
    pd.options.mode.chained_assignment = None  # default='warn'
    agg_snv_global_nonzero = agg_snv_global[cols].loc[:, agg_snv_global[cols].sum() > 0]
    cols = agg_snv_global_nonzero.columns.tolist()
    agg_snv_global_nonzero = pd.concat([agg_snv_global_nonzero, agg_snv_global[meta_cols]], axis=1)
    dynamic_data = []
    # return agg_snv_global_nonzero
    for patient, data in agg_snv_global_nonzero.groupby('draw_id'):
        p_draw_times = data.draw_times.unique().tolist()
        # p_draw_times.sort()
        if len(p_draw_times) >=3: #0 in p_draw_times and 1 in p_draw_times and 2 in p_draw_times:
            data_for_slope = data #[meta_cols + cols]
            for col in cols:
                slope, intercept, r, p, se = scipy.stats.linregress(data_for_slope.draw_month, data_for_slope[col])
                data_for_slope[f'{col}_slope'] = slope
                data[f'{col}_intercept'] = intercept
                # data[col] = data_for_slope[col]
                # data[:, f'{col}_r'] = r
                # data[:, f'{col}_p'] = p
            dynamic_data.append(data_for_slope)
    return pd.concat(dynamic_data, axis = 0)

def mutsig_global_to_df(mutsig_path_prefix, VCF_BASE_DIR) -> pd.DataFrame:
    input_vcf_paths2 = glob.glob(os.path.join(VCF_BASE_DIR, '*.vcf'))
    # input_vcf_paths2 = glob.glob(os.path.join(VCF_BASE_DIR, 'full_genome', '*.ann.filtered[1,3].vcf'))
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




