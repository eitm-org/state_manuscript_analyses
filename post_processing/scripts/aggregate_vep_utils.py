import os
import glob

import pandas as pd
import numpy as np
import warnings

from constants import refs_dir
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

LOCATION_INDEX = 1
FEATURE_TYPE_INDEX = 5
CONSEQUENCE_INDEX = 6
EXTRA_INDEX = 13

def get_vep_data_for_vcf(vcf_path):
    vep_data = {
        'chr': [],
        'position': [],
        'feature_type': [],
        # 'consequence': [],
        'consequence_impact': [],
        'distance_to_transcript': [],
        'sift_class': [],
        'sift_score': [],
        'polyphen_class': [],
        'polyphen_score': [],
        'am_class': [],
        'am_pathogenicity': [],
    }
    vcf_file_read = open(vcf_path, 'r')
    for line in vcf_file_read:
        if line.startswith('#'):
            continue
        line = line[1:-1].split('\t')
        location = line[LOCATION_INDEX]
        feature_type = line[FEATURE_TYPE_INDEX]
        consequence = line[CONSEQUENCE_INDEX]
        extra_field = line[EXTRA_INDEX]
        chr, position = location.split(':')
        vep_data['chr'].append(chr)
        vep_data['position'].append(int(position))
        vep_data['feature_type'].append(feature_type)
        # vep_data['consequence'].append(consequence)
        vep_data = parse_extra(extra_field, vep_data)
    vep_df = pd.DataFrame(vep_data)
    # return vep_df
    polyphen_class = vep_df.pop('polyphen_class')
    polyphen_class_ohe = pd.get_dummies(polyphen_class, prefix = 'polyphen_class')
    sift_class = vep_df.pop('sift_class')
    sift_class_ohe = pd.get_dummies(sift_class, prefix = 'sift_class')
    feature_type = vep_df.pop('feature_type')
    feature_type_ohe = pd.get_dummies(sift_class, prefix = 'feature_type')
    consequence_impact = vep_df.pop('consequence_impact')
    consequence_impact_ohe = pd.get_dummies(sift_class, prefix = 'consequence_impact')
    am_class = vep_df.pop('am_class')
    am_class_ohe = pd.get_dummies(am_class, prefix = 'am_class')
    return pd.concat([vep_df, polyphen_class_ohe, sift_class_ohe, feature_type_ohe, consequence_impact_ohe, am_class_ohe], axis = 1)

def parse_extra(extra_field, vep_data):
    extra = extra_field.split(';')
    extra_mapping = {e.split('=')[0]: e.split('=')[1] for e in extra}
    vep_data['consequence_impact'].append(extra_mapping.get('IMPACT'))
    vep_data['distance_to_transcript'].append(float(extra_mapping.get('DISTANCE') or np.nan))
    vep_data['am_class'].append(extra_mapping.get('am_class') or None)
    vep_data['am_pathogenicity'].append(float(extra_mapping.get('am_pathogenicity') or np.nan))
    if extra_mapping.get('SIFT'):
        sift_class = extra_mapping['SIFT'].split('(')[0]
        sift_score = float(extra_mapping['SIFT'].split('(')[-1][:-1])
    else:
        sift_class = None
        sift_score = np.nan
    if extra_mapping.get('PolyPhen'):
        polyphen_class = extra_mapping['PolyPhen'].split('(')[0]
        polyphen_score = float(extra_mapping['PolyPhen'].split('(')[-1][:-1])
    else:
        polyphen_class = None
        polyphen_score = np.nan
    vep_data['sift_class'].append(sift_class)
    vep_data['sift_score'].append(sift_score)
    vep_data['polyphen_class'].append(polyphen_class)
    vep_data['polyphen_score'].append(polyphen_score)
    return vep_data

def aggregate_vep_global(vep_df):
    vep_df_cp = vep_df.copy()
    vep_df_cp.pop('chr')
    vep_df_cp.pop('position')
    vep_df_agg = vep_df_cp.groupby('EIBS').sum()
    vep_df_agg  = vep_df_agg.add_suffix('_count')
    vep_df_cont_agg = vep_df_cp.groupby('EIBS').agg(
        {
            'distance_to_transcript': ['mean', 'min', 'max', 'std'],
            'sift_score': ['mean', 'min', 'max', 'std'],
            'am_pathogenicity': ['mean', 'min', 'max', 'std'],
            'polyphen_score': ['mean', 'min', 'max', 'std'],
        }
    )
    vep_df_cont_agg.columns = vep_df_cont_agg.columns.map('_'.join)
    return vep_df_agg.join(vep_df_cont_agg)

def aggregate_vep_chrom(vep_df):
    vep_df_cp = vep_df.copy()
    vep_df_cp.pop('position')
    vep_df_agg = vep_df_cp.groupby(['EIBS', 'chr']).sum()
    vep_df_agg  = vep_df_agg.add_suffix('_count')
    vep_df_cont_agg = vep_df_cp.groupby(['EIBS', 'chr']).agg(
        {
            'distance_to_transcript': ['mean', 'min', 'max', 'std'],
            'sift_score': ['mean', 'min', 'max', 'std'],
            'am_pathogenicity': ['mean', 'min', 'max', 'std'],
            'polyphen_score': ['mean', 'min', 'max', 'std'],
        }
    )
    vep_df_cont_agg.columns = vep_df_cont_agg.columns.map('_'.join)
    return vep_df_agg.join(vep_df_cont_agg)

def aggregate_vep_Mbps(vep_df):
    hg38_bed_path = os.path.join(refs_dir, 'GRCh38/GRCh38.primary_assembly.genome_X_chr.bed')
    hg38_bed = pd.read_csv(hg38_bed_path, sep='\t', header=None)
    hg38_bed.columns = ['chr', 'start', 'end']
    vep_df_cp = vep_df.copy()
    
    vep_df_cp['start'] = np.nan
    vep_df_cp['end'] = np.nan
    for _, (chr, start, end) in hg38_bed.iterrows():
        for i in range(start, end, int(1e6)):
            window_min = i
            window_max = min(end, int(i+3e6))
            vep_df_cp.loc[(vep_df['position'] > window_min) & (vep_df_cp['position'] <= window_max) & (vep_df_cp['chr'] == chr), 'start'] = window_min
            vep_df_cp.loc[(vep_df['position'] > window_min) & (vep_df_cp['position'] <= window_max) & (vep_df_cp['chr'] == chr), 'end'] = window_max       
    vep_df_cp.pop('position')
    vep_df_agg = vep_df_cp.groupby(['EIBS', 'chr', 'start', 'end']).sum()
    vep_df_cont_agg = vep_df_cp.groupby(['EIBS', 'chr', 'start', 'end']).agg(
        {
            'distance_to_transcript': ['mean', 'min', 'max', 'std'],
            'sift_score': ['mean', 'min', 'max', 'std'],
            'am_pathogenicity': ['mean', 'min', 'max', 'std'],
            'polyphen_score': ['mean', 'min', 'max', 'std'],
        }
    )
    vep_df_cont_agg.columns = vep_df_cont_agg.columns.map('_'.join)
    return vep_df_agg.join(vep_df_cont_agg)


def get_vep_df(vep_dir):
    vep_dfs = []
    for vcf_path in glob.glob(os.path.join(vep_dir, '*.vcf')):
        eibs = vcf_path.split('/')[-1].split('_vs_')[0]
        vep_df = get_vep_data_for_vcf(vcf_path)
        vep_df['EIBS'] = eibs
        vep_dfs.append(vep_df)
    return pd.concat(vep_dfs, axis = 0)

