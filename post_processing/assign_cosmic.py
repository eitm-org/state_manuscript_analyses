import os

import pandas as pd
from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh38', rsync=False, bash=True)

from STATE_analyses.post_processing.utils import scratch_dir

def process_sample_names(df, col):
    df['EIBS'] = df[col]
    parsed_name = df[col].str.split('_vs_', expand=True)
    df[col] = parsed_name[0]
    return df


def main():
    filternums = ['f1', 'f3']
    for filternum in filternums:
        vcf_path = os.path.join(scratch_dir, f'STATE_vcfs_{filternum}_region_filtered', 'full_genome')
        # TODO: remove input, output, logs dirs from vcf path
        Analyze.cosmic_fit(samples = vcf_path, output = f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/mut_sig_{filternum}', input_type = 'vcf', context_type="96",
                        genome_build="GRCh38", collapse_to_SBS96 = False)
        assignment_solution_activities = pd.read_csv(f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/mut_sig_{filternum}/Assignment_Solution/Activities/Assignment_Solution_Activities.txt', sep = '\t')
        assignment_solution_activities = process_sample_names(assignment_solution_activities, 'Samples') 
        assignment_solution_activities = assignment_solution_activities.drop(columns = ['EIBS'])
        assignment_solution_activities.rename(columns = {'Samples': 'EIBS'}, inplace=True)
        sbs_columns = [col for col in assignment_solution_activities.columns if col.startswith('SBS')]
        for col in sbs_columns:
            assignment_solution_activities = assignment_solution_activities.rename(columns={col: f'{col}_count'})
        assignment_solution_activities.to_csv(f'/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/mut_sig_assignment_solution_activities_{filternum}.csv', index = False)

    hg002_fp_vcfs = '/data/scratch/xchen/STATE_HG002_vcfs/fp'
    Analyze.cosmic_fit(samples = hg002_fp_vcfs, output = '/home/xchen@okta-oci.eitm.org/projects/STATE_analyses/data/mut_sig_fp', input_type = 'vcf', context_type="96",
                   genome_build="GRCh38", collapse_to_SBS96 = False)
if __name__ == "__main__":
    main()