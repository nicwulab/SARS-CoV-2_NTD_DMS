#!/usr/bin/python
import sys
import pandas as pd
import numpy as np

def count_to_freq(df, colname):
    df[colname+'_freq'] = (df[colname]+1)/(df[colname].sum()+len(df))
    return (df)

def exp_score_calculate(df, rep):
    exp_weight = df['Expr_bin0_'+rep+'_freq']*0.25 + df['Expr_bin1_'+rep+'_freq']*0.5 + \
                 df['Expr_bin2_'+rep+'_freq']*0.5  + df['Expr_bin3_'+rep+'_freq']*1
    exp_norm_factor = df['Expr_bin0_'+rep+'_freq'] + df['Expr_bin1_'+rep+'_freq'] + \
                      df['Expr_bin2_'+rep+'_freq'] + df['Expr_bin3_'+rep+'_freq']
    df['Exp_weight_'+rep] = exp_weight/exp_norm_factor
    w_summary = df.groupby('mut_class')['Exp_weight_'+rep].mean()
    w_summary = w_summary.reset_index()
    w_silent   = (float(w_summary.loc[w_summary['mut_class']=='silent']['Exp_weight_'+rep]))
    w_nonsense = (float(w_summary.loc[w_summary['mut_class']=='nonsense']['Exp_weight_'+rep]))
    df['Exp_score_'+rep] = (df['Exp_weight_'+rep]-w_nonsense)/(w_silent-w_nonsense)
    return (df)

def fusion_score_calculate(df, rep):
    df['Fus_ratio_'+rep] = np.log10(df['Fus_pos_'+rep+'_freq']/df['Fus_neg_'+rep+'_freq'])
    fus_summary = df.groupby('mut_class')['Fus_ratio_'+rep].mean()
    fus_summary = fus_summary.reset_index()
    fus_silent   = (float(fus_summary.loc[fus_summary['mut_class']=='silent']['Fus_ratio_'+rep]))
    fus_nonsense = (float(fus_summary.loc[fus_summary['mut_class']=='nonsense']['Fus_ratio_'+rep]))
    df['Fus_score_'+rep] = (df['Fus_ratio_'+rep]-fus_nonsense)/(fus_silent-fus_nonsense)
    return (df)

def wrapper(count_file):
    print ('reading: %s' % count_file)
    df = pd.read_csv(count_file, sep='\t')
    colnames = [colname for colname in df]
    for colname in colnames:
        if 'mut' not in colname and 'frag' not in colname:
            df = count_to_freq(df, colname)
    df = exp_score_calculate(df, 'rep1')
    df = exp_score_calculate(df, 'rep2')
    df = fusion_score_calculate(df, 'rep1')
    df = fusion_score_calculate(df, 'rep2')
    df['Exp_score'] = (df['Exp_score_rep1'] + df['Exp_score_rep2'])/2
    df['Fus_score'] = (df['Fus_score_rep1'] + df['Fus_score_rep2'])/2
    return (df)

def main():
    outfile = "result/NTD_DMS_scores.tsv"
    df_A = wrapper("result/NTD_DMS_count_aa_A.tsv")
    df_B = wrapper("result/NTD_DMS_count_aa_B.tsv")
    df   = pd.concat([df_A, df_B])
    print ('writing: %s' % outfile)
    df.to_csv(outfile, sep="\t", index=False)

if __name__ == "__main__":
    main()
