import os

import pandas as pd
from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh38', rsync=False, bash=True)

from constants import flat_results_dir

def process_sample_names(df, col):
    df['EIBS'] = df[col]
    parsed_name = df[col].str.split('_vs_', expand=True)
    df[col] = parsed_name[0]
    return df


def main():
    if not os.path.exists('cosmics_sig_data'):
        os.makedirs('cosmics_sig_data')
    filternums = ['f3']
    for filternum in filternums:
        vcf_path = os.path.join(flat_results_dir, f'STATE_vcfs_{filternum}_region_filtered', 'full_genome')
        Analyze.cosmic_fit(samples = vcf_path, output = f'cosmics_sig_data/mut_sig_{filternum}', input_type = 'vcf', context_type="96",
                        genome_build="GRCh38", collapse_to_SBS96 = False)
        assignment_solution_activities = pd.read_csv(f'cosmics_sig_data/mut_sig_{filternum}/Assignment_Solution/Activities/Assignment_Solution_Activities.txt', sep = '\t')
        assignment_solution_activities = process_sample_names(assignment_solution_activities, 'Samples') 
        assignment_solution_activities = assignment_solution_activities.drop(columns = ['EIBS'])
        assignment_solution_activities.rename(columns = {'Samples': 'EIBS'}, inplace=True)
        sbs_columns = [col for col in assignment_solution_activities.columns if col.startswith('SBS')]
        for col in sbs_columns:
            assignment_solution_activities = assignment_solution_activities.rename(columns={col: f'{col}_count'})
        assignment_solution_activities = assignment_solution_activities.set_index('EIBS')
        assignment_solution_activities['snv_count'] = assignment_solution_activities.sum(axis=1)
        assignment_solution_activities = assignment_solution_activities.reset_index()
        assignment_solution_activities.to_csv(f'cosmics_sig_data/mut_sig_assignment_solution_activities_{filternum}.csv', index = False)

    hg002_fp_vcfs = os.path.join(flat_results_dir,'STATE_HG002_vcfs/fp')
    Analyze.cosmic_fit(samples = hg002_fp_vcfs, output = 'cosmics_sig_data/mut_sig_hg002_fp', input_type = 'vcf', context_type="96",
                   genome_build="GRCh38", collapse_to_SBS96 = False)
if __name__ == "__main__":
    main()