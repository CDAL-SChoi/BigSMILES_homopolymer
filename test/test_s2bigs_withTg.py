from s2bigs import SMILES2BigSMILES
import time

start = time.time()

k = SMILES2BigSMILES()

#k.Data_Load_csv('with_Tg/Bicerano_bigsmiles.csv', SMILES='SMILES', encoding='cp949')
#k.s2bigs_converting('with_Tg/Bicerano_')

k.Data_Load_csv('with_Tg/JCIM_sup_bigsmiles.csv', SMILES='SMILES', encoding='cp949')
k.s2bigs_converting('with_Tg/JCIM_', move_parallel=0)

print(time.time()-start)