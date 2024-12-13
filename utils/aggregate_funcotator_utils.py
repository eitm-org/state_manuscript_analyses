import os
import glob
import random

from cyvcf2 import VCF, Writer
import pandas as pd
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

hg38_bed_path = '/data/xchen/refs/GRCh38/GRCh38.primary_assembly.genome_X_chr.bed'
hg38_bed = pd.read_csv(hg38_bed_path, sep='\t', header=None)
hg38_bed.columns = ['chr', 'start', 'end']



def get_funcotator_data_for_vcf(vcf_path):
    SBS_mapping = {
        'A>C': 'T>G',
        'A>G': 'T>C',
        'A>T': 'T>A',
        'C>A': 'C>A',
        'C>G': 'C>G',
        'C>T': 'C>T',
        'G>A': 'C>T',
        'G>C': 'C>G',
        'G>T': 'C>A',
        'T>A':  'T>A',
        'T>C': 'T>C',
        'T>G': 'T>G',
    }
    vcf = VCF(vcf_path)
    for h in vcf.header_iter():
        if h.info().get('ID') and h.info().get('ID') == 'FUNCOTATION': break
    func_desc = h.info()['Description']
    func_desc_list = func_desc.split(': ')[-1].split('|')
    func_dicts = []
    SBS = []
    for rec in vcf:
        rec_func = dict(zip(func_desc_list, rec.INFO.get('FUNCOTATION')[1:-1].split('|')))
        rec_func['CGC_Mutation_Type'] = rec_func['CGC_Mutation_Type'].replace('_', '').replace('%2C%20', ',')
        func_dicts.append(rec_func)
        SBS.append(SBS_mapping['>'.join([rec.REF, rec.ALT[0]])])
    func_df = pd.DataFrame(func_dicts)
    func_df.columns = func_df.columns.str.split('_', n=1, expand=True)
    func_df.replace('', 'NULL', inplace=True)
    func_df['SBS'] = SBS
    return func_df


def aggregate_funcotator_global(vcf_dir, debug=False):
    aggregate_func_data = defaultdict(list)
    for vcf_path in glob.glob(os.path.join(vcf_dir, '*.vcf')):
        eibs = vcf_path.split('/')[-1].split('_vs_')[0]
        if debug: print(eibs)
        func_df = get_funcotator_data_for_vcf(vcf_path)
        func_df_gencode = func_df['Gencode']
        func_df_cgc = func_df['CGC']
        func_df_clinvar = func_df['ClinVar']
        func_df_dbsnp = func_df['dbSNP']
        func_df_uniprot = func_df['Simple']
        func_df_sbs = func_df['SBS']
        aggregate_func_data['EIBS'].append(eibs)
        aggregate_sbs(func_df_sbs, aggregate_func_data)
        aggregate_uniprot(func_df_uniprot, aggregate_func_data)
        aggregate_dbsnp(func_df_dbsnp, aggregate_func_data)
        aggregate_cgc(func_df_cgc, aggregate_func_data)
        aggregate_gencode(func_df_gencode, aggregate_func_data)
    return pd.DataFrame(aggregate_func_data)

def aggregate_funcotator_chr(vcf_dir, debug=False):
    aggregate_func_data = defaultdict(list)
    for vcf_path in glob.glob(os.path.join(vcf_dir, '*.vcf')):
        eibs = vcf_path.split('/')[-1].split('_vs_')[0]
        if debug: print(eibs)
        func_df = get_funcotator_data_for_vcf(vcf_path)
        for chr, func_df_chr in func_df.groupby(('Gencode', '34_chromosome')):
            aggregate_func_data['chr'].append(chr)
            aggregate_func_data['EIBS'].append(eibs)
            func_df_gencode = func_df_chr['Gencode']
            func_df_cgc = func_df_chr['CGC']
            func_df_clinvar = func_df_chr['ClinVar']
            func_df_dbsnp = func_df_chr['dbSNP']
            func_df_uniprot = func_df_chr['Simple']
            func_df_sbs = func_df_chr['SBS']
            aggregate_sbs(func_df_sbs, aggregate_func_data)
            aggregate_uniprot(func_df_uniprot, aggregate_func_data)
            aggregate_dbsnp(func_df_dbsnp, aggregate_func_data)
            aggregate_cgc(func_df_cgc, aggregate_func_data)
            aggregate_gencode(func_df_gencode, aggregate_func_data)
    return pd.DataFrame(aggregate_func_data)


def aggregate_funcotator_Mpbs(vcf_dir, debug=True):
    aggregate_func_data = defaultdict(list)
    for vcf_path in glob.glob(os.path.join(vcf_dir, '*.vcf')):
        eibs = vcf_path.split('/')[-1].split('_vs_')[0]
        if debug:
            print(eibs)
        func_df = get_funcotator_data_for_vcf(vcf_path)
        func_df[('Gencode', '34_start')] = func_df['Gencode']['34_start'].astype(int)
        for chr, func_df_chr in func_df.groupby(('Gencode', '34_chromosome')):
            start, end = int(hg38_bed.loc[hg38_bed.chr == chr, 'start'].iloc[0]), int(hg38_bed.loc[hg38_bed.chr == chr, 'end'].iloc[0])
            for i in range(start, end, int(1e6)):
                window_min = i
                window_max = min(end, int(i+3e6))
                aggregate_func_data['chr'].append(chr)
                aggregate_func_data['EIBS'].append(eibs)
                aggregate_func_data['start'].append(window_min + 1)
                aggregate_func_data['end'].append(window_max)
                func_df_chr_window = func_df_chr.loc[(func_df_chr['Gencode']['34_start'] > window_min) & (func_df_chr['Gencode']['34_start'] <= window_max)]
                func_df_gencode = func_df_chr_window['Gencode']
                func_df_cgc = func_df_chr_window['CGC']
                func_df_clinvar = func_df_chr_window['ClinVar']
                func_df_dbsnp = func_df_chr_window['dbSNP']
                func_df_uniprot = func_df_chr_window['Simple']
                func_df_sbs = func_df_chr_window['SBS']
                aggregate_sbs(func_df_sbs, aggregate_func_data)
                aggregate_uniprot(func_df_uniprot, aggregate_func_data)
                aggregate_dbsnp(func_df_dbsnp, aggregate_func_data)
                aggregate_cgc(func_df_cgc, aggregate_func_data)
                aggregate_gencode(func_df_gencode, aggregate_func_data)
    return pd.DataFrame(aggregate_func_data)

def aggregate_sbs(func_df_sbs, aggregate_func_data):
    aggregate_func_data['sbs_ct_count'].append((func_df_sbs == 'C>T').sum())
    aggregate_func_data['sbs_cg_count'].append((func_df_sbs == 'C>G').sum())
    aggregate_func_data['sbs_ca_count'].append((func_df_sbs == 'C>A').sum())
    aggregate_func_data['sbs_ta_count'].append((func_df_sbs == 'T>A').sum())
    aggregate_func_data['sbs_tc_count'].append((func_df_sbs == 'T>C').sum())
    aggregate_func_data['sbs_tg_count'].append((func_df_sbs == 'T>G').sum())

def aggregate_uniprot(func_df_uniprot, aggregate_func_data):
    aggregate_func_data['uniprot_accession_exists_count'].append((func_df_uniprot['Uniprot_uniprot_accession'] != 'NULL').sum())

def aggregate_dbsnp(func_df_dbsnp, aggregate_func_data):
    aggregate_func_data['dbsnp_id_exists_count'].append((func_df_dbsnp['ID'] != 'NULL').sum())


def aggregate_cgc(func_df_cgc, aggregate_func_data):
    aggregate_func_data['cgc_entrez_id_count'].append((func_df_cgc['GeneID'] != 'NULL').sum())
    aggregate_func_data['cgc_entrez_unknown_count'].append((func_df_cgc['GeneID'] == 'NULL').sum())
    aggregate_func_data['cgc_somatic_in_gene_count'].append((func_df_cgc['Cancer_Somatic_Mut'] != 'NULL').sum())
    aggregate_func_data['cgc_germline_in_gene_count'].append((func_df_cgc['Cancer_Germline_Mut'] != 'NULL').sum())
    aggregate_func_data['cgc_missense_count'].append((func_df_cgc['Mutation_Type'].str.contains('Mis')).sum())
    aggregate_func_data['cgc_amplification_count'].append((func_df_cgc['Mutation_Type'].str.contains('A')).sum())
    aggregate_func_data['cgc_large_del_count'].append((func_df_cgc['Mutation_Type'].str.contains('D')).sum())
    aggregate_func_data['cgc_frameshift_count'].append((func_df_cgc['Mutation_Type'].str.contains('F')).sum())
    aggregate_func_data['cgc_nonsense_count'].append((func_df_cgc['Mutation_Type'].str.contains('N')).sum())
    aggregate_func_data['cgc_splicesite_count'].append((func_df_cgc['Mutation_Type'].str.contains('S')).sum())
    aggregate_func_data['cgc_translocation_count'].append((func_df_cgc['Mutation_Type'].str.contains('S')).sum())
    aggregate_func_data['cgc_unknown_count'].append((func_df_cgc['Mutation_Type'] == 'NULL').sum())


def aggregate_gencode(func_df_gencode, aggregate_func_data):
    # hugo
    aggregate_func_data['gencode_hugo_gene_count'].append((func_df_gencode['34_hugoSymbol'] != 'Unknown').sum())
    aggregate_func_data['gencode_hugo_unknown_count'].append((func_df_gencode['34_hugoSymbol'] == 'Unknown').sum())
    # variant class
    aggregate_func_data['gencode_variantclass_could_not_determine_count'].append((func_df_gencode['34_variantClassification'] == 'COULD_NOT_DETERMINE').sum())
    aggregate_func_data['gencode_variantclass_intron_count'].append((func_df_gencode['34_variantClassification'] == 'INTRON').sum())
    aggregate_func_data['gencode_variantclass_five_prime_utr_count'].append((func_df_gencode['34_variantClassification'] == 'FIVE_PRIME_UTR').sum())
    aggregate_func_data['gencode_variantclass_three_prime_utr_count'].append((func_df_gencode['34_variantClassification'] == 'THREE_PRIME_UTR').sum())
    aggregate_func_data['gencode_variantclass_igr_count'].append((func_df_gencode['34_variantClassification'] == 'IGR').sum())
    aggregate_func_data['gencode_variantclass_five_prime_flank_count'].append((func_df_gencode['34_variantClassification'] == 'FIVE_PRIME_FLANK').sum())
    aggregate_func_data['gencode_variantclass_three_prime_flank_count'].append((func_df_gencode['34_variantClassification'] == 'THREE_PRIME_FLANK').sum())
    aggregate_func_data['gencode_variantclass_missense_count'].append((func_df_gencode['34_variantClassification'] == 'MISSENSE').sum())
    aggregate_func_data['gencode_variantclass_nonsense_count'].append((func_df_gencode['34_variantClassification'] == 'NONSENSE').sum())
    aggregate_func_data['gencode_variantclass_nonstop_count'].append((func_df_gencode['34_variantClassification'] == 'NONSTOP').sum())
    aggregate_func_data['gencode_variantclass_silent_count'].append((func_df_gencode['34_variantClassification'] == 'SILENT').sum())
    aggregate_func_data['gencode_variantclass_splice_site_count'].append((func_df_gencode['34_variantClassification'] == 'SPLICE_SITE').sum())
    aggregate_func_data['gencode_variantclass_in_frame_del_count'].append((func_df_gencode['34_variantClassification'] == 'IN_FRAME_DEL').sum())
    aggregate_func_data['gencode_variantclass_in_frame_ins_count'].append((func_df_gencode['34_variantClassification'] == 'IN_FRAME_INS').sum())
    aggregate_func_data['gencode_variantclass_frame_shift_ins_count'].append((func_df_gencode['34_variantClassification'] == 'FRAME_SHIFT_INS').sum())
    aggregate_func_data['gencode_variantclass_frame_shift_del_count'].append((func_df_gencode['34_variantClassification'] == 'FRAME_SHIFT_DEL').sum())
    aggregate_func_data['gencode_variantclass_start_codon_snp_count'].append((func_df_gencode['34_variantClassification'] == 'START_CODON_SNP').sum())
    aggregate_func_data['gencode_variantclass_start_codon_ins_count'].append((func_df_gencode['34_variantClassification'] == 'START_CODON_INS').sum())
    aggregate_func_data['gencode_variantclass_start_codon_del_count'].append((func_df_gencode['34_variantClassification'] == 'START_CODON_DEL').sum())
    aggregate_func_data['gencode_variantclass_de_novo_start_in_frame_count'].append((func_df_gencode['34_variantClassification'] == 'DE_NOVO_START_IN_FRAME').sum())
    aggregate_func_data['gencode_variantclass_de_novo_start_out_frame_count'].append((func_df_gencode['34_variantClassification'] == 'DE_NOVO_START_OUT_FRAME').sum())
    aggregate_func_data['gencode_variantclass_rna_count'].append((func_df_gencode['34_variantClassification'] == 'RNA').sum())
    aggregate_func_data['gencode_variantclass_lincrna_count'].append((func_df_gencode['34_variantClassification'] == 'LINCRNA').sum())
    # transcript strand
    aggregate_func_data['gencode_transcriptstrand_plus_strand_count'].append((func_df_gencode['34_transcriptStrand'] == '+').sum())
    aggregate_func_data['gencode_transcriptstrand_minus_strand_count'].append((func_df_gencode['34_transcriptStrand'] == '-').sum())
    aggregate_func_data['gencode_transcriptstrand_null_strand_count'].append((func_df_gencode['34_transcriptStrand'] == 'NULL').sum())
    # gc
    aggregate_func_data['gencode_gc_min'].append(func_df_gencode['34_gcContent'].astype(float).min())
    aggregate_func_data['gencode_gc_avg'].append(func_df_gencode['34_gcContent'].astype(float).mean())
    aggregate_func_data['gencode_gc_max'].append(func_df_gencode['34_gcContent'].astype(float).max())
    # protein and cdna change
    aggregate_func_data['genecode_proteinchange_count'].append((func_df_gencode['34_proteinChange'] != 'NULL').sum())
    aggregate_func_data['genecode_cdnachange_count'].append((func_df_gencode['34_cDnaChange'] != 'NULL').sum())
    