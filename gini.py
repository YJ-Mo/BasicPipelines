import numpy as np
import pandas as pd
import re
from numba import jit
from scipy.stats import ks_2samp
import statsmodels.stats.multitest as multi

jit(nopython=False)
def get_gini(R_, nucleotide_, cut_AC=False, ratio_nan=0.9, ratio_nan_AC=0.8):
    ########
    if len(R_)==0:
        gini=np.nan
    else:
        R_nonull = R_[np.isnan(R_) == False]
        # R_nonull = R_[R_>0]
        ratio = len(R_nonull) / len(R_)
        if ratio <= ratio_nan:
            gini = np.nan
        else:
            if cut_AC:
                R_AC = R_[(nucleotide_ == b'A') | (nucleotide_ == b'C')]
                R_AC_nonull = R_AC[np.isnan(R_AC) == False]
                ratio_AC = len(R_AC_nonull) / len(R_AC)
                if (ratio_AC <= ratio_nan_AC) | len(R_AC) <= 1 | len(R_AC_nonull) <= 1:
                    gini = np.nan
                else:
                    sorted = np.sort(R_AC_nonull)
                    height, area = 0, 0
                    for i in range(0, len(sorted)):
                        height += sorted[i]
                        area += height - sorted[i] / 2.
                    fair_area = height * len(sorted) / 2.
                    if fair_area == 0:
                        gini = np.nan
                    else:
                        gini = (fair_area - area) / fair_area
            else:
                sorted = np.sort(R_nonull)
                height, area = 0, 0
                for i in range(0, len(sorted)):
                    height += sorted[i]
                    area += height - sorted[i] / 2.
                fair_area = height * len(sorted) / 2.
                if fair_area == 0:
                    gini = np.nan
                else:
                    gini = (fair_area - area) / fair_area
    return gini

@jit(nopython=False)
def get_p(R_1,R_2,ratio_nan=0.9):
    R_1_nonull = R_1[np.isnan(R_1) == False]
    ratio_1 = len(R_1_nonull) / len(R_1)
    R_2_nonull = R_2[np.isnan(R_2) == False]
    ratio_2 = len(R_2_nonull) / len(R_2)
    if (ratio_1 <= ratio_nan)|(ratio_2 <= ratio_nan):
        p = np.nan
    else:
        s,p=ks_2samp(R_1, R_2)
    return p

@jit(nopython=False)
def get_window(X, window=50, step=1, cut_AC=False, ratio_nan=0.9, ratio_nan_AC=0.8):
    #########读取数据#############
    name = X[0]

    nucleotide = np.array(list(X[1]))

    modified = np.array(X[2].split(',')).astype('float').astype('int')


    modified_depth = np.array(X[3].split(',')).astype('float').astype('int')

    #
    unmodified = np.array(X[4].split(',')).astype('float').astype('int')


    unmodified_depth = np.array(X[5].split(',')).astype('float').astype('int')
    Mutrs = modified / modified_depth
    Mutru = unmodified / unmodified_depth
    # R = Mutrs - Mutru
    R = np.array(X[6].split(',')).astype('float')
    R[R == -999] = np.nan
    n = len(nucleotide)

    #####滑动窗口计算change##################
    nucleotide_ = np.zeros([len(range(0, n - (window), step)), window], dtype=np.string_)
    R_ = np.zeros([len(range(0, n - (window), step)), window])
    gini = np.zeros([len(range(0, n - (window), step)), ])
    j = 0
    for i in range(0, n - (window), step):
        nucleotide_[j, :] = nucleotide[i:i + window]
        R_[j, :] = R[i:i + window]
        gini[j] = get_gini(R[i:i + window], nucleotide[i:i + window], cut_AC=cut_AC, ratio_nan=ratio_nan,
                           ratio_nan_AC=ratio_nan_AC)
        j = j + 1

    return nucleotide_, R_, gini

@jit(nopython=False)
def get_window_p(X1,X2, window=50, step=1,ratio_nan=0.9):
    #########读取数据#############
    name = X1[0]
    nucleotide = np.array(list(X1[1]))
    R_z = np.array(X1[6].split(',')).astype('float')
    R_z[R_z == -999] = np.nan
    R_f = np.array(X2[6].split(',')).astype('float')
    R_f[R_f == -999] = np.nan
    n = len(nucleotide)

    #####滑动窗口计算change##################
    nucleotide_ = np.zeros([len(range(0, n - (window), step)), window], dtype=np.string_)
    R_z_ = np.zeros([len(range(0, n - (window), step)), window])
    R_f_ = np.zeros([len(range(0, n - (window), step)), window])
    p = np.zeros([len(range(0, n - (window), step)), ])
    j = 0
    for i in range(0, n - (window), step):
        nucleotide_[j, :] = nucleotide[i:i + window]
        R_z_[j, :] = R_z[i:i + window]
        R_f_[j, :] = R_f[i:i + window]
        p[j] = get_p(R_z[i:i + window], R_f[i:i + window],ratio_nan=ratio_nan,)
        j = j + 1
    p_=p.copy()
    if len(p[~np.isnan(p)])>0:
        a,p_bh,b,c =multi.multipletests(p[~np.isnan(p)],method='fdr_bh')
        p_[~np.isnan(p_)]=p_bh
    return p_,p


@jit(nopython=False)
def calculate_delta_gini(R_1, R_2, gini_1, gini_2, ratio_nan=0.9):
    '''Calculates Standard RMSD on two vectors of numbers of the same length'''
    # Check to see the vectors are of equal length.
    #计算delta_gini,(1)当R_1,R_2长度不同时，输出nan值;(2)当滑动窗口非空值比率<90%，gini_index输出空值。(3)除(1)(2)情况外，输出delta_gini。
    if len(R_1) != len(R_2):
        return np.nan
    else:
        R = R_1 - R_2
        if len(R[np.isnan(R) == False]) / len(R) <= ratio_nan:
            return np.nan
        else:
            delta_gini = gini_1 - gini_2
    return delta_gini

@jit(nopython=False)
def get_all(X,cut_AC=False, ratio_nan=0, ratio_nan_AC=0.8):
    #########读取数据#############
    name = X[0]
    nucleotide = np.array(list(X[1]))
    R = np.array(X[6].split(',')).astype('float')
    R[R == -999] = np.nan
    n = len(nucleotide)
    gini= get_gini(R, nucleotide, cut_AC=cut_AC, ratio_nan=ratio_nan,ratio_nan_AC=ratio_nan_AC)
    return gini

@jit(nopython=False)
def get_se(X,location ,start,end,strand,cut_AC=False, ratio_nan=0, ratio_nan_AC=0.8):
    #########读取数据#############
    name = X[0]
    nucleotide = np.array(list(X[1]))
    R = np.array(X[6].split(',')).astype('float')
    R[R == -999] = np.nan
    n = len(nucleotide)
    if strand=='+':
        R_=R[location<start]
        nucleotide_=nucleotide[location<start]
        gini_5= get_gini(R_, nucleotide_, cut_AC=cut_AC, ratio_nan=ratio_nan,ratio_nan_AC=ratio_nan_AC)
        R_ = R[location > end]
        nucleotide_ = nucleotide[location > end]
        gini_3 = get_gini(R_, nucleotide_, cut_AC=cut_AC, ratio_nan=ratio_nan, ratio_nan_AC=ratio_nan_AC)
    else:
        R_ = R[location < start]
        nucleotide_ = nucleotide[location < start]
        gini_3 = get_gini(R_, nucleotide_, cut_AC=cut_AC, ratio_nan=ratio_nan, ratio_nan_AC=ratio_nan_AC)
        R_ = R[location > end]
        nucleotide_ = nucleotide[location > end]
        gini_5 = get_gini(R_, nucleotide_, cut_AC=cut_AC, ratio_nan=ratio_nan, ratio_nan_AC=ratio_nan_AC)
    R_ = R[(location<= end)&(location>=start)]
    nucleotide_ = nucleotide[(location<= end)&(location>=start)]
    gini_cds = get_gini(R_, nucleotide_, cut_AC=cut_AC, ratio_nan=ratio_nan, ratio_nan_AC=ratio_nan_AC)
    return gini_3,gini_cds,gini_5


@jit(nopython=False)
def get_statistics(data_z, data_f,data_gff_,location_list,strand):
    data_z = np.array(data_z).reshape([7, ])
    data_f = np.array(data_f).reshape([7, ])
    nucleotide_z, R_z, gini_z = get_window(data_z)
    nucleotide_f, R_f, gini_f = get_window(data_f)
    ###uv+###
    gini_z_=gini_z[pd.notnull(gini_z)]
    if len(gini_z_)==0:
        gini_max_z=np.nan
    else:
        gini_max_z=gini_z[pd.notnull(gini_z)].max()
    gini_all_z=get_all(data_z)
    if len(data_gff_.loc[data_gff_['location']=='CDS','strat'])==0:
        gini_3_z=np.nan
        gini_5_z=np.nan
        gini_cds_z=np.nan
    else:
        CDS_strat=list(data_gff_.loc[data_gff_['location']=='CDS','strat'])[0]
        CDS_end = list(data_gff_.loc[data_gff_['location']=='CDS','end'])[0]
        gini_3_z, gini_cds_z, gini_5_z = get_se(data_z,location_list,CDS_strat,CDS_end,strand)
    ###uv-###
    gini_f_ = gini_f[pd.notnull(gini_f)]
    if len(gini_f_) == 0:
        gini_max_f = np.nan
    else:
        gini_max_f = gini_f[pd.notnull(gini_f)].max()
    gini_all_f = get_all(data_f)
    if len(data_gff_.loc[data_gff_['location'] == 'CDS', 'strat']) == 0:
        gini_3_f = np.nan
        gini_5_f = np.nan
        gini_cds_f = np.nan
    else:
        CDS_strat = list(data_gff_.loc[data_gff_['location'] == 'CDS', 'strat'])[0]
        CDS_end = list(data_gff_.loc[data_gff_['location'] == 'CDS', 'end'])[0]
        gini_3_f, gini_cds_f, gini_5_f = get_se(data_f, location_list, CDS_strat, CDS_end, strand)

    if len(R_z) != len(R_f):
        print("error")
    else:
        delta_gini = np.zeros([len(R_z), ])
        for i in range(len(R_z)):
            delta_gini[i] = calculate_delta_gini(R_z[i, :], R_f[i, :], gini_z[i], gini_f[i])
    return delta_gini,gini_max_z,gini_all_z,gini_3_z, gini_cds_z, gini_5_z,gini_max_f,gini_all_f,gini_3_f, gini_cds_f, gini_5_f


def main_sum_gini(data_uv_z, data_uv_f, data_location, data_gff, transcript_id_list):
    statistics_sum = pd.DataFrame(columns={'transcript_id', 'num', 'num_0.1','delta_max', 'delta_min'})
    statistics_z = pd.DataFrame(columns={'transcript_id', 'gini_all', 'gini_max', 'gini_3UTR', 'gini_5UTR', 'gini_CDS'})
    statistics_f = pd.DataFrame(columns={'transcript_id', 'gini_all', 'gini_max', 'gini_3UTR', 'gini_5UTR', 'gini_CDS'})
    i = 0
    for transcript in transcript_id_list:
        data_z = data_uv_z.loc[data_uv_z['transcript_id'] == transcript, :]
        data_f = data_uv_f.loc[data_uv_f['transcript_id'] == transcript, :]
        data_location_ = data_location.loc[data_location['transcript_id'] == transcript,:]
        data_gff_ = data_gff.loc[data_gff['transcript_id'] == transcript, :]
        location_list=np.array(list(data_location_['location'])[0].split(',')).astype('int')
        chr = list(data_location_['chr'])[0]
        strand=list(data_location_['strand'])[0]
        delta_gini, gini_max_z, gini_all_z, gini_3_z, gini_cds_z, gini_5_z, gini_max_f, gini_all_f, gini_3_f, gini_cds_f, gini_5_f= get_statistics(data_z,data_f, data_gff_, location_list, strand)

        statistics_sum.loc[i, 'transcript_id'] = transcript
        num = sum(np.abs(delta_gini) >=0.2)
        statistics_sum.loc[i, 'num'] = num
        num_1 = sum(np.abs(delta_gini) >=0.1)
        statistics_sum.loc[i, 'num_0.1'] = num_1
        statistics_sum.loc[i, 'delta_gini_list'] = ','.join(list(delta_gini.astype('str')))

        if len(delta_gini[np.isnan(delta_gini) == False]) > 0:
            statistics_sum.loc[i, 'delta_max'] = np.max(delta_gini[np.isnan(delta_gini) == False])
            statistics_sum.loc[i, 'delta_min'] = np.min(delta_gini[np.isnan(delta_gini) == False])
        else:
            statistics_sum.loc[i, 'delta_max'] = np.nan
            statistics_sum.loc[i, 'delta_min'] = np.nan

        statistics_z.loc[i,'transcript_id']=transcript
        statistics_z.loc[i,'gini_all']=gini_all_z
        statistics_z.loc[i, 'gini_max'] = gini_max_z
        statistics_z.loc[i, 'gini_3UTR'] = gini_3_z
        statistics_z.loc[i, 'gini_5UTR'] = gini_5_z
        statistics_z.loc[i, 'gini_CDS'] = gini_cds_z

        statistics_f.loc[i,'transcript_id']=transcript
        statistics_f.loc[i,'gini_all']=gini_all_f
        statistics_f.loc[i, 'gini_max'] = gini_max_f
        statistics_f.loc[i, 'gini_3UTR'] = gini_3_f
        statistics_f.loc[i, 'gini_5UTR'] = gini_5_f
        statistics_f.loc[i, 'gini_CDS'] = gini_cds_f
        i = i + 1
        if i % 100 == 0:
            print(i)
    return statistics_sum,statistics_z,statistics_f

if __name__ == '__main__':
    col_uv_f = 'your_dir/final.modified_umodified/col_nouv/'
    col_uv_z = 'your_dir/final.modified_umodified/col_uv/'


###col###
    result = 'your_output_dir/'
    uv_z = col_uv_z
    uv_f = col_uv_f

    data_uv_z = pd.read_csv(uv_z + '/final.modified_unmodified_new', sep='\t')
    data_uv_z = data_uv_z[['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                         'Untreated_mutations', 'Untreated_effective_depth', 'reactivity']]
    data_uv_z = data_uv_z.drop_duplicates()

    data_uv_f = pd.read_csv(uv_f + '/final.modified_unmodified_new', sep='\t')
    data_uv_f = data_uv_f[['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                         'Untreated_mutations', 'Untreated_effective_depth', 'reactivity']]
    data_uv_f = data_uv_f.drop_duplicates()
    transcript_uv_z = pd.read_csv(uv_z + '/cutoff.hit.group', sep='\t')
    transcript_uv_z.columns = ['cut_off', 'transcript_id', 'modified_depth_median',
                               'unmodified_depth_median', 'modified_depth_sum', 'unmodified_depth_sum', 'hit_level']

    transcript_uv_f = pd.read_csv(uv_f + '/cutoff.hit.group', sep='\t')
    transcript_uv_f.columns = ['cut_off', 'transcript_id', 'modified_depth_median',
                               'unmodified_depth_median', 'modified_depth_sum', 'unmodified_depth_sum', 'hit_level']

    uv_z_transcript = list(transcript_uv_z.loc[(transcript_uv_z['modified_depth_median'] > 100) & (
            transcript_uv_z['hit_level'] > 0), 'transcript_id'])

    uv_f_transcript = list(transcript_uv_f.loc[(transcript_uv_f['modified_depth_median'] > 100) & (
            transcript_uv_f['hit_level'] > 0), 'transcript_id'])
    transcript_all = list(set(uv_z_transcript) & set(uv_f_transcript))
    print(len(transcript_all))
    gff_path = 'your_dir/Ath_genes.gff'
    data_gff = pd.read_csv(gff_path, sep='\t')
    data_location = pd.read_csv(
        'your_dir/transcript_exon_location.csv',
        sep='\t')
    statistics_sum, statistics_z, statistics_f = main_sum_gini(data_uv_z, data_uv_f, data_location, data_gff, transcript_all)
    pd.DataFrame(statistics_sum).to_csv(result + '/gini_summary_50_1.csv', sep='\t', index=False, header=True)
    pd.DataFrame(statistics_z).to_csv(result + '/gini_summary_UV+_50_1.csv', sep='\t', index=False, header=True)
    pd.DataFrame(statistics_f).to_csv(result + '/gini_summary_UV-_50_1.csv', sep='\t', index=False, header=True)
