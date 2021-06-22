## calculating transcript numbers at varying hit levels
import numpy as np
import pandas as pd
import re
from numpy import median
from numba import jit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


# path
exon_gtf_path='your_dir/ath_exons.gtf'
col_uv_f_path='your_dir/final.modified_umodified/col_nouv/'
col_uv_z_path='your_dir/final.modified_umodified/col_uv/'

# extract info from references
gtf_data=pd.read_csv(exon_gtf_path,sep='\t',header=None)
gtf_data_new=pd.DataFrame(columns={'transcript_id','gene_type','gene_name','chr'})
gtf_data_new['transcript_id']=gtf_data.iloc[:,8].apply(lambda x:x.split(';')[1].split('"')[1])
gtf_data_new['gene_type'] = gtf_data.iloc[:,8].apply(lambda x:re.findall('gene_biotype ".*?"',x)[0].split('"')[1])
gtf_data_new['gene_name'] = gtf_data.iloc[:,8].apply(lambda x:re.findall('gene_name ".*?"',x)[0].split('"')[1] if 'gene_name' in x else np.nan)
gtf_data_new['chr'] = gtf_data.iloc[:,0].apply(lambda x: 6 if x=='Mt' else 7 if x=='Pt' else x ).astype('int')
gtf_data_new = gtf_data_new.drop_duplicates()
gtf_data_new.index = range(len(gtf_data_new))

# merge extracted ref info with output hit data
hit_level_col_uv_f = pd.read_csv(col_uv_f_path+'/final.modified_unmodified.hit',sep='\t',header=None)
hit_level_col_uv_f.columns =['transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit']
hit_level_col_uv_f = pd.merge(hit_level_col_uv_f,gtf_data_new,on='transcript_id',how='left')
hit_level_col_uv_f['Type'] = 'WT_UV-'

hit_level_col_uv_z = pd.read_csv(col_uv_z_path+'/final.modified_unmodified.hit',sep='\t',header=None)
hit_level_col_uv_z.columns =['transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit']
hit_level_col_uv_z = pd.merge(hit_level_col_uv_z,gtf_data_new,on='transcript_id',how='left')
hit_level_col_uv_z['Type'] = 'WT_UV+'

# calculate transcript numbers at different hit level threshold
data1=pd.DataFrame(columns={'Type','Number of transcripts','hit'})
for num in [-1000,0,1,2,5,10,15]:
    WT_UV_f = len(hit_level_col_uv_f.loc[(hit_level_col_uv_f['hit'] > num) & (hit_level_col_uv_f['modified.median'] > 100), :])
    WT_UV_z = len(hit_level_col_uv_z.loc[(hit_level_col_uv_z['hit'] > num) & (hit_level_col_uv_z['modified.median'] > 100)])

    data2 = pd.DataFrame(columns={'Type', 'Number of transcripts', 'hit'})
    data2['Type'] = ['WT_UV-', 'WT_UV+']
    data2['Number of transcripts'] = [WT_UV_f, WT_UV_z]

    if num==-1000:
        data2['hit'] ='ALL transcripts'
    else:
        data2['hit']='Hit level>'+str(num)
    data1 = pd.concat([data1,data2])

# plot hit level distributions
plt.switch_backend('agg')
plt.figure(figsize=(8, 6))
sns.set(style="ticks", context="talk")
sns.axes_style({'font.family': ['sans-serif'],'font.sans-serif': ['Arial']})
g = sns.barplot(y='Number of transcripts',x='hit',hue='Type',data=data1,palette=["#3498DB", "#1F618D"])

sns.despine()
font1 = {'family' : 'Arial','weight' : 'roman','size': 22}
plt.xticks(rotation=60)

plt.legend(fontsize='small')
plt.tight_layout()
plt.savefig('your_dir/transcript_num_hit.png')
plt.close()
