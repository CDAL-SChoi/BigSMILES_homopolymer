from bigs2s import BigSMILES2SMILES
import time

start = time.time()

k = BigSMILES2SMILES()

k.Data_Load_csv('without_Tg/polyBERT_len85_0.csv', BigSMILES='0')

k.Converting('bigs2s_')

print(time.time()-start)
