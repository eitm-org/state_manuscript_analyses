import glob
import os

import pandas as pd
pd.set_option('mode.chained_assignment', None)

from constants import flat_results_dir

class RangeDict(dict):
    def __getitem__(self, item):
        if not isinstance(item, range): # or xrange in Python 2
            for key in self:
                if item in key:
                    return self[key]
            raise KeyError(item)
        else:
            return super().__getitem__(item) # or super(RangeDict, self) for Python 2

def get_noisy_coverage_bedfile(mosdepth_file):
    df = pd.read_csv(mosdepth_file, compression='gzip', sep='\t', header=None, names = ['chr', 'start', 'end', 'depth'])
    df = df.loc[(df.chr.str.startswith('chr')) & (df.chr != 'chrM')]
    df['chr_num'] = pd.to_numeric(df.chr.str[3:].replace('X', '23').replace('Y', '24'))
    df.sort_values(['chr_num','start'], inplace=True)
    window = df[['chr', 'depth']].groupby('chr').rolling(200).agg(['min', 'max', 'std'])
    window = window.droplevel(0)
    df['rolling_min'] = window[('depth', 'min')]
    df['rolling_max'] = window[('depth', 'max')]
    df['rolling_std'] = window[('depth', 'std')]

    df_autosome = df.loc[~ (df.chr.str.endswith('X') | df.chr.str.endswith('Y'))]
    df_sexsome = df.loc[df.chr.str.endswith('X') | df.chr.str.endswith('Y')]
    
    autosome_depth_thresh_min, autosome_depth_thresh_max = df_autosome.agg({'depth': [lambda x: x.quantile(0.05), lambda x: x.quantile(0.995)]}).depth.tolist()
    autosome_std_thresh = df_autosome.agg({'rolling_std': lambda x: x.quantile(0.995)}).tolist()[0]

    sexsome_depth_thresh_min, sexsome_depth_thresh_max = df_sexsome.agg({'depth': [lambda x: x.quantile(0.05), lambda x: x.quantile(0.995)]}).depth.tolist()
    sexsome_std_thresh = df_sexsome.agg({'rolling_std': lambda x: x.quantile(0.995)}).tolist()[0]

    df_autosome['depth_thresh_min'] = [autosome_depth_thresh_min] * df_autosome.shape[0]
    df_autosome['depth_thresh_max'] = [autosome_depth_thresh_max] * df_autosome.shape[0]
    df_autosome['std_thresh'] = [autosome_std_thresh] * df_autosome.shape[0]
    df_sexsome['depth_thresh_min'] = [sexsome_depth_thresh_min] * df_sexsome.shape[0] 
    df_sexsome['depth_thresh_max'] = [sexsome_depth_thresh_max] * df_sexsome.shape[0]
    df_sexsome['std_thresh'] = [sexsome_std_thresh] * df_sexsome.shape[0]

    df = pd.concat([df_autosome, df_sexsome], axis=0)
    df['window_filter'] = ( ((df.rolling_min > df.depth_thresh_min) &\
                                        (df.rolling_max < df.depth_thresh_max))  &\
                                        (df.rolling_std <  df.std_thresh) 
                                     )
    df.drop(columns='chr_num', inplace=True)
    return df.loc[~df.window_filter]

def main():
    mos_base_path = os.path.join(flat_results_dir, 'STATE_mosdepth')
    mosdepth_paths = glob.glob(os.path.join(mos_base_path, 'EIBS*'))
    noisy_cov_base_path = os.path.join(flat_results_dir, 'noisy_coverage_regions')
    if not os.path.exists(noisy_cov_base_path): 
        os.makedirs(noisy_cov_base_path)
    for mosdepth_file in mosdepth_paths:
        eibs = mosdepth_file.split('/')[-1].split('.')[0]
        noisy_cov_file = os.path.join(noisy_cov_base_path, f'{eibs}_noisy_cov.bed')
        if os.path.exists(noisy_cov_file): continue
        print(f'Noisy regions identified for {eibs}')
        noisy_cov = get_noisy_coverage_bedfile(mosdepth_file)
        noisy_cov.to_csv(noisy_cov_file, index=False, sep='\t')

if __name__ == "__main__":
    main()
