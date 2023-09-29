from bigs2s import BigSMILES2SMILES
import time

start = time.time()

k = BigSMILES2SMILES()

#k.Data_Load_csv('with_Tg/Bicerano_bigsmiles.csv', BigSMILES='BigSMILES', encoding='cp949')
#k.bigs2s_converting('with_Tg/bigs2s_journal')

k.Data_Load_csv('with_Tg/JCIM_sup_bigsmiles.csv', BigSMILES='BigSMILES', encoding='cp949')
k.bigs2s_converting('with_Tg/bigs2s_JCIM')

print(time.time()-start)