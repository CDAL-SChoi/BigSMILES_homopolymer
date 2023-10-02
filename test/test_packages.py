from BigSMILES_homopolymer import SMILES2BigSMILES as s2bigs

test=s2bigs()
smileslist=['*CCCCO*', '*N[Si](*)(C)C', '*CC=C(C*)CCCC', '*CC(*)OCCCCCCCC']
result=[]
for i in range(len(smileslist)):
  result.append(test.Converting_single(SMILES=smileslist[i]))

print(result)


from BigSMILES_homopolymer import BigSMILES2SMILES as bigs2s

test=bigs2s()
bigsmileslist=['{<CCCCO>}','{<N[Si](C)(C)>}','{$CC=C(CCCC)C$}','{$CC(OCCCCCCCC)$}']

result=[]
for i in range(len(bigsmileslist)):
  result.append(test.Converting_single(BigSMILES=bigsmileslist[i]))

print(result)
