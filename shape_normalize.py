import numpy as np
import pandas as pd
import re
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import median
from numba import jit

#########转录本注释和数据读入#################
#### path #####
exon_gtf_path='/data/TA_QUIZ_RNA_regulation/data/ATH/GTF/shape_map/ath_exons.gtf'
col_uv_f_path='/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_nouv/'
col_uv_z_path='/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_uv/'

### 读取注释文件，在exon_gtf文件的第8列中提取'transcript_id','gene_type','gene_name','chr'等信息
gtf_data=pd.read_csv(exon_gtf_path,sep='\t',header=None)
gtf_data_new=pd.DataFrame(columns={'transcript_id','gene_type','gene_name','chr'})
gtf_data_new['transcript_id']=gtf_data.iloc[:,8].apply(lambda x:x.split(';')[1].split('"')[1])
gtf_data_new['gene_type'] = gtf_data.iloc[:,8].apply(lambda x:re.findall('gene_biotype ".*?"',x)[0].split('"')[1])
gtf_data_new['gene_name'] = gtf_data.iloc[:,8].apply(lambda x:re.findall('gene_name ".*?"',x)[0].split('"')[1] if 'gene_name' in x else np.nan)
gtf_data_new['chr'] = gtf_data.iloc[:,0].apply(lambda x: 6 if x=='Mt' else 7 if x=='Pt' else x ).astype('int')
gtf_data_new = gtf_data_new.drop_duplicates()
gtf_data_new.index = range(len(gtf_data_new))

#  提取col_nouv、col_uv中'cutoff.hit.group'的'group','transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit'列，并将ath_exons参考基因组中的'gene_type'，'gene_name','chr'等信息合并到原数据中

hit_level_col_uv_f = pd.read_csv(col_uv_f_path+'/cutoff.hit.group',sep='\t')
hit_level_col_uv_f.columns =['group','transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit']
hit_level_col_uv_f = pd.merge(hit_level_col_uv_f,gtf_data_new,on='transcript_id',how='left')
hit_level_col_uv_f['spe'] = 'WT_UV-'

hit_level_col_uv_z = pd.read_csv(col_uv_z_path+'/cutoff.hit.group',sep='\t')
hit_level_col_uv_z.columns =['group','transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit']
hit_level_col_uv_z = pd.merge(hit_level_col_uv_z,gtf_data_new,on='transcript_id',how='left')
hit_level_col_uv_z['spe'] = 'WT_UV+'


###从col_uv_f、col_uv_z中筛选出gene_type属于['lncRNA','rRNA','tRNA']，hit_level的值>10 && modified.median>5000的转录本，并输出它的转录本ID###
col_uv_f_id = hit_level_col_uv_f.loc[(hit_level_col_uv_f.gene_type.isin(['lncRNA','rRNA','tRNA']))&(hit_level_col_uv_f['modified.median']>5000)&(hit_level_col_uv_f['hit']>10),'transcript_id']
col_uv_z_id = hit_level_col_uv_z.loc[(hit_level_col_uv_z.gene_type.isin(['lncRNA','rRNA','tRNA']))&(hit_level_col_uv_z['modified.median']>5000)&(hit_level_col_uv_z['hit']>10),'transcript_id']

print(col_uv_f_id)
print(col_uv_z_id)
print(set(col_uv_f_id)&set(col_uv_z_id))
