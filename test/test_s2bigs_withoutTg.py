from s2bigs import SMILES2BigSMILES
import time

start = time.time()

k = SMILES2BigSMILES()

k.Data_Load_txt('dataset/generated_polymer_smiles_train.txt', N=10000000)
### This data is accessible at https://zenodo.org/record/7766806 (last checked 23.09.30)

k.Converting('polyBERT_len85_', move_parallel=0)

print(time.time()-start)
