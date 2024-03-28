import bigsmiles as bs

from BigSMILES_homopolymer import version_converter

import pandas as pd
import numpy as np

import logging
import os

logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)    ### to avoid too many logging messages


conv = version_converter()

count = 0
tot = 0

path = ['with_Tg/Bicerano_bigsmiles.csv', 'with_Tg/JCIM_sup_bigsmiles.csv']
## path = ['without_Tg/polyBERT_len85_0.csv']


for dataname in path:
    
    data = pd.read_csv(dataname, encoding='CP949')
    print(dataname)
    
    tot = tot + len(data)
    
    data = data['BigSMILES'].values   ## if with_Tg
    ## data = data['0'].values        ## if without_Tg
    
    
    for i in range(len(data)):
        tmp = conv.zero2one(data[i])
        try:
            bs.BigSMILES(tmp)
        except:
            count = count+1
            #print(tmp)      ## to see which data arise error
            continue

    
print('Total length')
print(tot)

print('Valid counts')
print(count)
