import numpy as np
import pandas as pd
import re
from numpy import median
from numba import jit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import math




## input ref exons/transcripts and data folder
exon_gtf_path='your_dir/ath_exons.gtf'
col_uv_f_path='your_dir/final.modified_umodified/col_nouv/'
col_uv_z_path='your_dir/final.modified_umodified/col_uv/'

## extract ref info from ref
gtf_data=pd.read_csv(exon_gtf_path,sep='\t',header=None)
gtf_data_new=pd.DataFrame(columns={'transcript_id','gene_type','gene_name','chr'})
gtf_data_new['transcript_id']=gtf_data.iloc[:,8].apply(lambda x:x.split(';')[1].split('"')[1])
gtf_data_new['gene_type'] = gtf_data.iloc[:,8].apply(lambda x:re.findall('gene_biotype ".*?"',x)[0].split('"')[1])
gtf_data_new['gene_name'] = gtf_data.iloc[:,8].apply(lambda x:re.findall('gene_name ".*?"',x)[0].split('"')[1] if 'gene_name' in x else np.nan)
gtf_data_new['chr'] = gtf_data.iloc[:,0].apply(lambda x: 6 if x=='Mt' else 7 if x=='Pt' else x ).astype('int')
gtf_data_new = gtf_data_new.drop_duplicates()
gtf_data_new.index = range(len(gtf_data_new))

##  extract data info from data
hit_level_col_uv_f = pd.read_csv(col_uv_f_path+'/cutoff.hit.group',sep='\t')
hit_level_col_uv_f.columns =['group','transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit']
hit_level_col_uv_f = pd.merge(hit_level_col_uv_f,gtf_data_new,on='transcript_id',how='left')
hit_level_col_uv_f['spe'] = 'WT_UV-'

hit_level_col_uv_z = pd.read_csv(col_uv_z_path+'/cutoff.hit.group',sep='\t')
hit_level_col_uv_z.columns =['group','transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit']
hit_level_col_uv_z = pd.merge(hit_level_col_uv_z,gtf_data_new,on='transcript_id',how='left')
hit_level_col_uv_z['spe'] = 'WT_UV+'

## select gene type, hit level and modified.median
col_uv_f_id = hit_level_col_uv_f.loc[(hit_level_col_uv_f.gene_type.isin(['lncRNA','rRNA','tRNA']))&(hit_level_col_uv_f['modified.median']>5000)&(hit_level_col_uv_f['hit']>10),'transcript_id']
col_uv_z_id = hit_level_col_uv_z.loc[(hit_level_col_uv_z.gene_type.isin(['lncRNA','rRNA','tRNA']))&(hit_level_col_uv_z['modified.median']>5000)&(hit_level_col_uv_z['hit']>10),'transcript_id']

print(col_uv_f_id)
print(col_uv_z_id)
print(set(col_uv_f_id)&set(col_uv_z_id))

@jit(nopython=False)
def calc_quartile(x, q, qtype=7):
    # 计算给定的四分位数，x为所需计算的四分位数的列表，q代表是第25%百分位数还是75%百分位数。qtype代表所选用的方法。
    # source: http://adorio-research.org/wordpress/?p=125
    # x = array, q = quartile (in % as a decimal)
    y = np.copy(x)
    n = len(y)
    abcd = [(0, 0, 1, 0),  # inverse empirical distrib.function., R type 1
            (0.5, 0, 1, 0),  # similar to type 1, averaged, R type 2
            (0.5, 0, 0, 0),  # nearest order statistic,(SAS) R type 3
            (0, 0, 0, 1),  # California linear interpolation, R type 4
            (0.5, 0, 0, 1),  # hydrologists method, R type 5
            (0, 1, 0, 1),  # mean-based estimate(Weibull method), (SPSS,Minitab), type 6
            (1, -1, 0, 1),  # mode-based method,(S, S-Plus), R type 7
            (1.0 / 3, 1.0 / 3, 0, 1),  # median-unbiased ,  R type 8
            (3 / 8.0, 0.25, 0, 1)  # normal-unbiased, R type 9.
            ]
    a, b, c, d = abcd[qtype - 1]
    g, j = math.modf(a + (n + b) * q - 1)
    #第1四分位数的计算公式为 （n+3)/4,这里为(n-1)/4;第3四分位数的计算公式为(3n+1)/4,这里为(3n-3)/4。是由于数组中的元素是从x[0]开始，因此需要-1.
    if j < 0:
        return x[0]
    elif j >= n:
        return x[n - 1]
    j = int(math.floor(j))
    if g == 0:
        return x[j]
    else:
        return y[j] + (y[j + 1] - y[j]) * (c + d * g)

# @jit(nopython=False)
def find_boxplot_factor(array):
    x, o, a = [], [], 0
    # Following deprecated line is behavior that normalization and
    # structure modeling were optimized with, but this behavior
    # is probably not ideal. For RNAs with regions of poor sequencing
    # depth, treating those regions as unreactive could bias the
    # normalized reactivities. This is especially important for
    # larger RNAs.
    # x = np.fromiter((n if not isnan(n) else 0 for n in array))
    x = array[np.where(np.isfinite(array))]
    # 除去数组中的非有限数的值。
    x = x[np.nonzero(x)]
    # 输出不为0的元素的下标
    if x.shape[0] < 10:
        norm_factor = np.nan
    else:
        x.sort()
        ten_pct = len(x) // 10
        five_pct = len(x) // 20
        # calculate the interquartile range *1.5
        q_limit = 1.5 * abs(calc_quartile(x, 0.25) - calc_quartile(x, 0.75))
        ten_limit = x[x.shape[0] - 1 - ten_pct]
        five_limit = x[x.shape[0] - 1 - five_pct]
        # choose the cutoff that eliminates the fewest points
        limit = max(q_limit, ten_limit)
        if len(x) < 100:
            limit = max(q_limit, five_limit)
        # make new list without the outliers
        for i in range(len(x)):
            if x[i] < limit:
                o.append(x[i])
        # avg next ten percent
        try:
            for i in range(-ten_pct, 0):
                a = o[i] + a
            norm_factor = a / ten_pct
        except IndexError:
            norm_factor = np.nan
    return norm_factor

@jit(nopython=False)
def find_boxplot_factor_new(array):
    ####the mean of the 92-98th percentile reactivities#####
    # x, o, a = [], [], []
    x = array[np.where(np.isfinite(array))]
    x = x[np.nonzero(x)]
    if x.shape[0] < 10:
        norm_factor = np.nan
    else:

        # calculate the interquartile range *1.5
        o = x[(x>=np.quantile(x,0.92))&(x<=np.quantile(x,0.98))]
        norm_factor = np.mean(o)
    return norm_factor


@jit(nopython=False)
def filter(X,Mutru_cut_off=0.02,read_depth_cutoff=100,R_window=10,window_mutru_cutoff=0.03,window_mutru_num=3,window_mutrs_cutoff=0.1,window_mutrs_num=3):
    name = X[0]
    nucleotide = X[1]
    modified = np.array(X[2].split(',')).astype('float').astype('int')
    modified_depth = np.array(X[3].split(',')).astype('float').astype('int')
    unmodified = np.array(X[4].split(',')).astype('float').astype('int')
    unmodified_depth = np.array(X[5].split(',')).astype('float').astype('int')

    # print(nucleotide)
    Mutrs = modified/modified_depth
    Mutru = unmodified/unmodified_depth
    R = Mutrs - Mutru
    n = len(modified)
    ##### 校正reactivity###########
    R[Mutru>Mutru_cut_off] = np.nan
    R[(modified_depth <= read_depth_cutoff)|(unmodified_depth)<=read_depth_cutoff] = np.nan
    R[R<0]= np.nan
    for i in range(len(R)-R_window+1):
        data_Mutru = Mutru[i:i+R_window-1]
        data_Mutrs = Mutrs[i:i+R_window-1]
        if (len(data_Mutru[data_Mutru>window_mutru_cutoff])>=window_mutru_num)|(len(data_Mutrs[data_Mutrs>window_mutrs_cutoff])>=window_mutrs_num):
            R[i:i + R_window-1] = np.nan
    R_new = ",".join(list(R.astype('str'))).replace('nan','-999')
    X = X.append(pd.Series(R_new))
    return X
    #
@jit(nopython=False)
def get_fatcor(data):
    R_all = np.array(data.iloc[0,6].split(',')).astype('float')
    for i in range(1,len(data)):
        R_ = np.array(data.iloc[i,6].split(',')).astype('float')
        R_all = np.hstack((R_all,R_))
    R_all[R_all==-999]=np.nan
    factor = find_boxplot_factor_new(R_all)
    return factor

@jit(nopython=False)
def normalization_all(X,factor):
    #########读取数据#############
    R = np.array(X[6].split(',')).astype('float')
    R_nan = len(R[R == -999])/len(R)
    R_0 = len(R[R == 0])/len(R)
    R[R == -999] = np.nan
    R_new = R/factor
    R_new = ",".join(list(R_new.astype('str'))).replace('nan','-999')
    X[6] =R_new
    return X,R_nan,R_0
def main(data,transcript_list=[]):

    data_filter = pd.DataFrame(np.zeros([len(data),7]))
    data_filter.columns = ['transcript_id','Nucleotide' ,'Modified_mutations', 'Modified_effective_depth',
                             'Untreated_mutations', 'Untreated_effective_depth', 'reactivity']
    data_new = pd.DataFrame(np.zeros([len(data), 7]))
    data_new.columns = ['transcript_id','Nucleotide' ,'Modified_mutations', 'Modified_effective_depth',
                           'Untreated_mutations', 'Untreated_effective_depth', 'reactivity']
    print(len(data))
    for i in range(len(data)):
        if i %1000 == 0:
            print(i)
        X=filter(data.iloc[i,:])
        data_filter.iloc[i,:] = list(X)
    if len(transcript_list)> 0:
        factor = get_fatcor(data_filter.loc[data_filter.transcript_id.isin(transcript_list),:])
    else:
        factor = get_fatcor(data_filter)
    for i in range(len(data_filter)):
        if i %1000 == 0:
            print(i)
            # print data_
        data_new.iloc[i,:],R_nan,R_0 = normalization_all(data_filter.iloc[i,:],factor)
        # data_new.loc[i, 'ratio_nan']=R_nan
        # data_new.loc[i, 'ratio_0'] = R_0
    data_new = data_new[['transcript_id','Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                           'Untreated_mutations', 'Untreated_effective_depth', 'reactivity']]
    return data_new

def loctaion(X):
    X_1 = X.split('(')[1].split(')')[0].split(',')
    loctaion_list = []
    for i in range(len(X_1)):
        loctaion_list_ = [i for i in range(int(X_1[i].split(':')[0]),int(X_1[i].split(':')[1])+1)]
        loctaion_list.extend(loctaion_list_)
    loctaion_list =np.array(loctaion_list)
    loctaion_list = ",".join(list(loctaion_list.astype('str')))

    return loctaion_list



def mapping(exon_data):
    data_location = pd.DataFrame(np.zeros([len(exon_data),4]))
    data_location.columns = ['transcript_id','chr','strand','location']

    for i in range(len(exon_data)):
        if i%1000 ==0:
            print(i)
        data_location.loc[i, 'transcript_id'] = exon_data.loc[i, 'transcript_id']
        data_location.loc[i,'chr'] = exon_data.loc[i,'chr'].split('.')[0]
        data_location.loc[i, 'strand'] = exon_data.loc[i, 'chr'].split('.')[1]
        data_location.loc[i, 'location'] = loctaion(exon_data.loc[i, 'start_end'])
    return data_location

if __name__ == '__main__':
    col_uv_f = 'your_dir/final.modified_umodified/col_nouv/'
    col_uv_z = 'your_dir/final.modified_umodified/col_uv/'
    col_uv_f2 = 'your_output_dir/final.modified_umodified/col_nouv/'
    col_uv_z2 = 'your_output_dir/final.modified_umodified/col_uv/'
    path=[col_uv_f,col_uv_z]
    path2=[col_uv_f2,col_uv_z2]
    for i in range(len(path)):
        transcript_list = [ 'AT3G41768.1', 'AT3G06355.1', 'ATMG01390.1' ]
        data_uv_z = pd.read_csv(path[i] + '/final.modified_unmodified', sep='\t')
        data_uv_z.columns = ['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth','Untreated_mutations', 'Untreated_effective_depth']
        data_uv_z_new=main(data_uv_z,transcript_list)
        print(data_uv_z_new)
        data_uv_z_new.to_csv(path2[i]+'/final.modified_unmodified_new',sep='\t',header=True,index=False)
