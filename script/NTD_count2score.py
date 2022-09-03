#!/usr/bin/python
import sys
import pandas as pd
import numpy as np

def count_to_freq(df, colname):
    df[colname+'_freq'] = (df[colname]+1)/(df[colname].sum()+len(df))
    return (df)
def count_to_total_freq(df):
    total_count_rep_1 = df['Expr_bin0_rep1'].sum() + df['Expr_bin1_rep1'].sum() + df['Expr_bin2_rep1'].sum()  + df['Expr_bin3_rep1'].sum()
    each_mutant_count_rep_1 = df['Expr_bin0_rep1'] + df['Expr_bin1_rep1'] + df['Expr_bin2_rep1']  + df['Expr_bin3_rep1']
    df['totalfreq_rep_1'] = each_mutant_count_rep_1/total_count_rep_1
    total_count_rep_2 = df['Expr_bin0_rep2'].sum() + df['Expr_bin1_rep2'].sum() + df['Expr_bin2_rep2'].sum()  + df['Expr_bin3_rep2'].sum()
    each_mutant_count_rep_2 = df['Expr_bin0_rep2'] + df['Expr_bin1_rep2'] + df['Expr_bin2_rep2']  + df['Expr_bin3_rep2']
    df['totalfreq_rep_2'] = each_mutant_count_rep_2/total_count_rep_2
    

    df['avg_total_freq'] = (df['totalfreq_rep_1'] + df['totalfreq_rep_2']) / 2
    return(df)

def exp_score_calculate(df, rep, freq_cutoff):
    print (rep)
    exp_weight = df['Expr_bin0_'+rep+'_freq']*0.25 + df['Expr_bin1_'+rep+'_freq']*0.5 + \
                 df['Expr_bin2_'+rep+'_freq']*0.75  + df['Expr_bin3_'+rep+'_freq']*1
    exp_norm_factor = df['Expr_bin0_'+rep+'_freq'] + df['Expr_bin1_'+rep+'_freq'] + \
                      df['Expr_bin2_'+rep+'_freq'] + df['Expr_bin3_'+rep+'_freq']
    df['Exp_weight_'+rep] = exp_weight/exp_norm_factor
    df_high_freq = df[df['avg_total_freq'] >= freq_cutoff]
    print(len(df_high_freq))
    w_summary = df_high_freq.groupby('mut_class')['Exp_weight_'+rep].mean()
    w_summary = w_summary.reset_index()
    print (w_summary)
    w_silent   = (float(w_summary.loc[w_summary['mut_class']=='silent']['Exp_weight_'+rep]))
    w_nonsense = (float(w_summary.loc[w_summary['mut_class']=='nonsense']['Exp_weight_'+rep]))
    df['Exp_score_'+rep] = (df['Exp_weight_'+rep]-w_nonsense)/(w_silent-w_nonsense)
    print ('w_nonsense', w_nonsense)
    print ('w_silent', w_silent)
    return (df)

def biomodality_analysis(df):
    df['sum_freq_rep1'] = df['Expr_bin0_rep1_freq'] + df['Expr_bin1_rep1_freq'] + df['Expr_bin2_rep1_freq']  + df['Expr_bin3_rep1_freq']
    #df['bin0_bin3_ratio_rep1'] = (df['Expr_bin0_rep1'] + df['Expr_bin3_rep1'])/each_mutant_count_rep_1
    df['bin0_ratio_rep1'] = df['Expr_bin0_rep1_freq']/df['sum_freq_rep1']
    df['bin3_ratio_rep1'] = df['Expr_bin3_rep1_freq']/df['sum_freq_rep1']

    df['sum_freq_rep2'] = df['Expr_bin0_rep2_freq'] + df['Expr_bin1_rep2_freq'] + df['Expr_bin2_rep2_freq']  + df['Expr_bin3_rep2_freq']
    #df['counts_across_bins_rep_2'] = each_mutant_count_rep_2
    #df['bin0_bin3_ratio_rep2'] = (df['Expr_bin0_rep2'] + df['Expr_bin3_rep2'])/each_mutant_count_rep_2
    df['bin0_ratio_rep2'] = df['Expr_bin0_rep2_freq']/df['sum_freq_rep2']
    df['bin3_ratio_rep2'] = df['Expr_bin3_rep2_freq']/df['sum_freq_rep2'] #using freq instead of counts for bimodality analysis
    return (df)

def fusion_score_calculate(df, rep, freq_cutoff):
    df['Normalized_fus_ratio_'+rep] = np.log10(df['Fus_pos_'+rep+'_freq']/df['Fus_neg_'+rep+'_freq']/df['avg_total_freq'])
    df_high_freq = df[df['Input_freq'] >= freq_cutoff]
    fus_summary = df_high_freq.groupby('mut_class')['Normalized_fus_ratio_'+rep].mean()
    fus_summary = fus_summary.reset_index()
    fus_silent   = (float(fus_summary.loc[fus_summary['mut_class']=='silent']['Normalized_fus_ratio_'+rep]))
    fus_nonsense = (float(fus_summary.loc[fus_summary['mut_class']=='nonsense']['Normalized_fus_ratio_'+rep]))
    df['Fus_score_'+rep] = (df['Normalized_fus_ratio_'+rep]-fus_nonsense)/(fus_silent-fus_nonsense)
    return (df)

def wrapper(count_file, freq_cutoff):
    print ('reading: %s' % count_file)
    df = pd.read_csv(count_file, sep='\t')
    colnames = [colname for colname in df]
    for colname in colnames:
        if 'mut' not in colname and 'frag' not in colname:
            df = count_to_freq(df, colname)
    df = count_to_total_freq(df)
    df = exp_score_calculate(df, 'rep1', freq_cutoff)
    df = exp_score_calculate(df, 'rep2', freq_cutoff)
    df['Exp_score'] = (df['Exp_score_rep1'] + df['Exp_score_rep2'])/2 
    df = fusion_score_calculate(df, 'rep1', freq_cutoff)
    df = fusion_score_calculate(df, 'rep2', freq_cutoff)
    df['Fus_score'] = (df['Fus_score_rep1'] + df['Fus_score_rep2'])/2
    df = biomodality_analysis(df)
    return (df)

def main():
    freq_cutoff = 0.000075
    outfile_1 = "result/NTD_DMS_scores.tsv"
    outfile_2 = "result/NTD_DMS_scores_by_resi.tsv"
    outfile_3 = "result/NTD_DMS_scores_by_resi_and_replicates.tsv"
    df_A = wrapper("result/NTD_DMS_count_aa_A.tsv", freq_cutoff)
    df_B = wrapper("result/NTD_DMS_count_aa_B.tsv", freq_cutoff)
    df   = pd.concat([df_A, df_B])
    print ('writing: %s' % outfile_1)
    df.to_csv(outfile_1, sep="\t", index=False)
    df['resi'] = df['mut'].str[0:-1]
    all_resi   = df[df['mut_class'] != 'WT'].groupby('resi').size().reset_index()
    df_by_resi = df
    df_by_resi = df_by_resi[df_by_resi['avg_total_freq'] >= freq_cutoff]
    df_by_resi = df_by_resi[df_by_resi['mut_class'] != 'WT']
    df_by_resi = df_by_resi[df_by_resi['mut_class'] != 'silent']
    df_by_resi = df_by_resi[df_by_resi['mut_class'] != 'nonsense']
    df_by_resi_mean  = df_by_resi.groupby('resi')['Exp_score'].mean().reset_index(name='mean_exp_score')
    df_by_resi_count = df_by_resi.groupby('resi').size().reset_index(name='count')
    df_by_resi = pd.merge(df_by_resi_mean, df_by_resi_count, on='resi', how='outer')
    df_by_resi = pd.merge(all_resi, df_by_resi, on='resi', how='outer')
    df_by_resi = df_by_resi.sort_values(by='resi', key=lambda x:x.str[1::].astype(int))
    df_by_resi['pos'] = df_by_resi['resi'].str[1::].astype(int)
    df_by_resi = df_by_resi[['resi','pos','count','mean_exp_score']]
    print ('writing: %s' % outfile_2)
    df_by_resi.to_csv(outfile_2, sep="\t", index=False)

    df_by_resi_and_replicates = df
    df_by_resi_and_replicates = df_by_resi_and_replicates[df_by_resi_and_replicates['avg_total_freq'] >= freq_cutoff]
    df_by_resi_and_replicates = df_by_resi_and_replicates[df_by_resi_and_replicates['mut_class'] != 'WT']
    df_by_resi_and_replicates = df_by_resi_and_replicates[df_by_resi_and_replicates['mut_class'] != 'silent']
    df_by_resi_and_replicates = df_by_resi_and_replicates[df_by_resi_and_replicates['mut_class'] != 'nonsense']
    df_by_resi_mean_rep1  = df_by_resi_and_replicates.groupby('resi')['Exp_score_rep1'].mean().reset_index(name='mean_exp_score_rep1')
    df_by_resi_mean_rep2  = df_by_resi_and_replicates.groupby('resi')['Exp_score_rep2'].mean().reset_index(name='mean_exp_score_rep2')
    df_by_resi_and_replicates = pd.merge(df_by_resi_mean_rep1, df_by_resi_mean_rep2, on='resi', how='outer')
    df_by_resi_and_replicates = pd.merge(df_by_resi_and_replicates, df_by_resi_count, on='resi', how='outer')
    df_by_resi_and_replicates = pd.merge(all_resi, df_by_resi_and_replicates, on='resi', how='outer')
    df_by_resi_and_replicates = df_by_resi_and_replicates.sort_values(by='resi', key=lambda x:x.str[1::].astype(int))
    df_by_resi_and_replicates['pos'] = df_by_resi_and_replicates['resi'].str[1::].astype(int)
    df_by_resi_and_replicates = df_by_resi_and_replicates[['resi','pos','count','mean_exp_score_rep1','mean_exp_score_rep2']]
    print ('writing: %s' % outfile_3)   
    df_by_resi_and_replicates.to_csv(outfile_3, sep="\t", index=False)

if __name__ == "__main__":
    main()