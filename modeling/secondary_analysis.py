'''
Assumes the data files are in a subdirectory called 'data/'
Assumes the evaluation-metric files are in a subdirectory called 'metrics/'
'''


import os
import glob
import itertools

import matplotlib
import matplotlib.pylab as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.stats import spearmanr


data_files = [
    'intent_affect_lagged_021623.csv',
    'urge_affect_lagged_021623.csv',
]


eval_metrics_cols = [
    'mae',
    'rmse',
    'r2',
]


comp_metrics_cols = [
    'percent_high',
    'percent_zero',    
    'marginal_entropy',
    'marginal_mean',
    'avg_time_diff_between_observations',
    'var_time_diff_between_observations',
    'num_observations',    
    'percent_small_jumps',
    'percent_big_jumps',
    'percent_high_to_entropy_ratio',
]


comp_metric_names = {
    'age': 'Age',
    'race_white': 'Race',
    'gender': 'Gender', 
    
    'percent_zero': 'Percent Zero',
    'percent_small_jumps': 'Percent\nSmall Jumps',
    'percent_big_jumps': 'Percent\nBig Jumps',
    'marginal_entropy': 'Marginal Entropy',
    'marginal_mean': 'Marginal Mean',
    'num_observations': 'Number of\nObservations',
    'avg_time_diff_between_observations': 'Avg Time Between\nObservations',
    'var_time_diff_between_observations': 'Var Time Between\nObservations',
    'r2': 'r2',
    'percent_high': 'Percent High',
    'percent_high_to_entropy_ratio': 'Percent High /\nMarginal Entropy',    
}


OUTPUT_DIR = 'output'


def csv_fname_to_title(fname):
    if 'ar1' in fname:
        model_name = 'Baseline Model (Simple Autoregressive Model)'
    elif 'gp_affect' in fname:
        model_name = 'Gaussian Process Model (with Affect Features)'
    elif 'gp' in fname:
        model_name = 'Simple Gaussian Process Model'
    else:
        model_name = 'Baseline Model (Elastic Net with Affect Features)'
        
    if 'intent' in fname:
        target_name = 'Predicting Suicidal Intent'
    elif 'urge' in fname:
        target_name = 'Predicting Suicidal Urges'
    else:
        raise Exception('Unknown target')
        
    return f'{model_name} {target_name}'


def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    demographics = pd.read_csv(os.path.join('data', 'demogs_included.csv'))
    
    for data_fname in data_files:
        outcome = data_fname.split('_')[0]
        df = pd.read_csv(os.path.join('data', data_fname))

        print('*' * 79)
        print(f'* OUTCOME = {outcome}:')
        print('*' * 79)
        print('')
        
        metric_files = glob.glob(f'metrics/*{outcome}*.csv')
        for metric_fname in metric_files:

            print(f'Results file: {metric_fname}')
            
            eval_metrics = pd.read_csv(metric_fname)
            assert(len(eval_metrics['ppt_id'].unique()) == len(eval_metrics))
            
            def compute_metrics(row):
                subdf = df.loc[df['ppt_id'] == row['ppt_id']]
                curr = subdf['sitb_si_' + outcome]
                prev = subdf[outcome + '_lead']
                diff = np.abs(curr - prev)
                
                row['residual_magnitude_variance'] = np.var(diff[~np.isnan(diff)])
                row['percent_zero'] = (curr == 0.0).mean() 
                row['percent_high'] = (curr >= 5.0).mean()
                row['percent_small_jumps'] = (diff <= 3.0).mean()
                row['percent_big_jumps'] = (diff > 3.0).mean()

                marginal_freq = curr.value_counts()
                for i in range(11):
                    if i not in marginal_freq:
                        marginal_freq.at[i] = 0
                    marginal_freq.at[i] = marginal_freq[i] + 1
                    
                marginal_cat = (marginal_freq / float(marginal_freq.sum()))
                row['marginal_entropy'] = -(marginal_cat * np.log(marginal_cat)).sum()
                
                row['marginal_mean'] = np.mean(curr)

                row['percent_high_to_entropy_ratio'] = row['percent_high'] / row['marginal_entropy']
                
                row['num_observations'] = len(curr)

                time_in_study = np.sort(subdf['time_in_study'].values)
                row['avg_time_diff_between_observations'] = np.mean(
                    (time_in_study[1:] - time_in_study[:-1]) ** 2.0
                )

                row['var_time_diff_between_observations'] = np.var(
                    (time_in_study[1:] - time_in_study[:-1]) ** 2.0
                )

                subdem = demographics.loc[demographics['ppt_id'] == row['ppt_id']]
                if len(subdem) == 0:
                    return row
                
                assert(len(subdem) == 1)
                dem = subdem.iloc[0]

                row['age'] = dem['Age']
                row['race_white'] = float(dem['Race_White'])
                row['gender'] = float(dem['GenderID'] == 'Male')

                eval_m = eval_metrics[eval_metrics['ppt_id'] == row['ppt_id']]
                assert(len(eval_m) == 1)

                row['r2'] = eval_m.iloc[0]['r2']
                
                return row

            comp_metrics = eval_metrics.apply(compute_metrics, axis='columns')

            for cm in comp_metrics_cols:
                print('Top for: {}'.format(comp_metric_names[cm]))
                subdf = comp_metrics.loc[comp_metrics[cm].nlargest(10).index]
                subdf = subdf[['ppt_id', cm]]
                print(subdf)
            
            results = pd.DataFrame()
            for em, cm in itertools.product(eval_metrics_cols, comp_metrics_cols):
                corr = spearmanr(comp_metrics[[em, cm]], nan_policy='omit')

                asterisk = ''
                if corr.pvalue < 0.001:
                    asterisk = '\star\star\star'
                elif corr.pvalue < 0.01:
                    asterisk = '\star\star'
                elif corr.pvalue < 0.05:
                    asterisk = '\star'

                if asterisk == '':
                    annotation = f'{corr.statistic:.2f}'
                else:
                    annotation = r'${' + f'{corr.statistic:.2f}' + r'}^{' + asterisk + r'}$'
                
                results = results.append({
                    'Evaluation Metric': em,
                    'Time Series Feature': comp_metric_names[cm],
                    'Correlation': corr.statistic,
                    'P-Value': corr.pvalue,
                    'Annotation': annotation,
                }, ignore_index=True)
            
            hm = results.pivot('Evaluation Metric', 'Time Series Feature', 'Correlation')
            am = results.pivot('Evaluation Metric', 'Time Series Feature', 'Annotation')

            b_size = 1.1
            fig, ax = plt.subplots(
                1, 1, figsize=(b_size * (len(comp_metrics_cols) + 1),
                               1.3 + b_size * len(eval_metrics_cols)),
            )
            
            sns.heatmap(hm, annot=am, ax=ax, fmt='', linewidth=0.5, vmin=-1.0, vmax=1.0,
                        cmap=sns.color_palette('coolwarm', as_cmap=True))

            ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
            plt.tight_layout()
            plt.title(csv_fname_to_title(metric_fname))
            plt.savefig(f'{OUTPUT_DIR}/{os.path.basename(metric_fname).split(".")[0]}.png', dpi=200)
            plt.close()
            
            print('\n')


if __name__ == '__main__':
    main()

