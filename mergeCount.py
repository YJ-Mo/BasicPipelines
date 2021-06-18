import numpy as np
import pandas as pd
import re
import os

sample=['sample_1','sample_2','sample_3']
path = 'your_dir/5.read_counts/result/'
data_list=[]
for i in range(len(sample)):
    print(sample[i])
    temp=pd.read_csv(path+sample[i]+'.all.txt',sep='\t',index_col=0)
    data_list.append(temp)
matrix_C = data_list[0]
for i in range(len(data_list)-1):
    matrix_C = matrix_C.join(data_list[i+1],how="outer")
matrix_C = matrix_C.sort_index(axis=1)

name=sample[0][0:3:1] #3 should be length of sample name
matrix_C.to_csv('your_dir/6.mergeCount/'+name+'_count_all.txt',sep='\t',index=True,header=True)
