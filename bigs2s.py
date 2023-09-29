import pandas as pd
import numpy as np

import csv
from openpyxl import load_workbook

import re
from collections import Counter

from rdkit import Chem


###         Sunho Choi (schoi_@korea.ac.kr)
###         Korea University, School of Electrical Engineering
###         23.09.30.
###
###         Article:
###         An automated BigSMILES conversion workflow and dataset for homopolymeric macromolecules
###



class BigSMILES2SMILES:
    def __init__(self):
        self.SMILES_data = None
        self.num = 0

        self.BigSMILES_data = None
        
        
    def Data_Load_xlsx(self, filename, data_only=True, BigSMILES='BigSMILES'):
        load_wb = load_workbook(filename, data_only)
        load_ws = load_wb.active
        data = pd.DataFrame([[i.value for i in j] for j in load_ws.rows])
        header = data.iloc[0]
        data = data[1:]
        data.rename(columns=header, inplace=True)
        data.dropna(subset=[BigSMILES], inplace=True)
        data.reset_index(drop=True, inplace=True)

        self.BigSMILES_data = data[BigSMILES]

        self.num = len(self.BigSMILES_data)
        
        
    def Data_Load_csv(self, filename, BigSMILES='BigSMILES', encoding='utf-8'):
        ## TXT format files can also be imported
        data = pd.read_csv(filename, header=0, encoding=encoding)
        data.dropna(subset=[BigSMILES], inplace=True)
        data.reset_index(drop=True, inplace=True)

        self.BigSMILES_data = data[BigSMILES]

        self.num = len(self.BigSMILES_data)

        
    def bigs2s_converting(self, resultfile='result'):
        self.SMILES_data = []
        kk = 0
        
        for i in range(self.num):
            
            print('process : ', i, ' / ', self.num)
            
            BigSMILES = self.BigSMILES_data[i]
            
            if '>,<' in BigSMILES:
                ind_x = BigSMILES.rfind(',')
                BigSMILES = BigSMILES[:ind_x]
                
            for j in ['<', '>', '{', '}', '\$', ',', ' ']:
                BigSMILES = re.sub(j,'',BigSMILES)

            BigSMILES = '*' + BigSMILES + '*'

            BigSMILES = Chem.MolToSmiles(Chem.MolFromSmiles(BigSMILES))
            
            self.SMILES_data.append(BigSMILES)
            
            if i - kk == 100000:
                aa = pd.DataFrame(self.SMILES_data)
                aa.to_csv(resultfile+str(kk)+'.csv')
                kk = i
                self.SMILES_data = []
                aa = None
                
        if kk != i:
            aa = pd.DataFrame(self.SMILES_data)
            aa.to_csv(resultfile+str(kk)+'.csv')
            kk = i
            self.SMILES_data = []
            aa = None
            
        

