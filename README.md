<img src="https://capsule-render.vercel.app/api?type=egg&color=99CCFF&height=50&section=header"/>


# Automated BigSMILES conversion workflow and dataset for homopolymeric macromolecules

This repository contains the code of generating BigSMILES from the SMILES representations of homopolymers, and the opposite process as well.
The dataset is available from the [Figshare](https://figshare.com/s/813ca7794bd9257e9843).

**Abstract:**
Simplified molecular-input line-entry system (SMILES) has been utilized in a variety of artificial intelligence analyses owing to its capability to represent chemical structures in line notation. However, its ease of representation is limited, which has led to the proposal of BigSMILES as an alternative method suitable for the representation of macromolecules. Nevertheless, the studies on BigSMILES remains limited due to its preprocessing requirements. Thus, this study proposes a conversion workflow of BigSMILES, focusing on its automated generation from SMILES representations of homopolymers. BigSMILES representations for 4,927,181 records are provided, therby enabling its immediate uses for various research and development applications. Our study presents detailed descriptions on a validation process to ensure the accuracy, interchangeability, and robustness of the conversion. Additionally, a systematic overview of codes and functions utilized is present, emphasizing their relevance in the context of BigSMILES generation. This advancement is anticipated to significantly aid researchers, facilitating further studies in BigSMILES representation, including potential applications in deep learning and further extension to complex structures such as copolymers.

## Install BigSMILES_homopolymer package
```py
pip install git+https://github.com/CDAL-SChoi/BigSMILES_homopolymer.git
```

## Example code for SMILES to BigSMILES conversion of single SMILES

```py
from BigSMILES_homopolymer import SMILES2BigSMILES as s2bigs

test=s2bigs()
smileslist=['*CCCCO*', '*N[Si](*)(C)C', '*CC=C(C*)CCCC', '*CC(*)OCCCCCCCC']
result=[]
for i in range(len(smileslist)):
  result.append(test.Converting_single(SMILES=smileslist[i]))

print(result)
```

## Example code for BigSMILES to SMILES conversion of single BigSMILES

```py
from BigSMILES_homopolymer import BigSMILES2SMILES as bigs2s

test=bigs2s()
bigsmileslist=['{<CCCCO>}','{<N[Si](C)(C)>}','{$CC=C(CCCC)C$}','{$CC(OCCCCCCCC)$}']

result=[]
for i in range(len(bigsmileslist)):
  result.append(test.Converting_single(BigSMILES=bigsmileslist[i]))

print(result)
```

## The workflow of the s2bigs.py code
![figure 3](https://github.com/CDAL-SChoi/BigSMILES_homopolymer/assets/50295574/ec6b6b03-1387-459f-b823-ba429bd91ebd)



<img src="https://capsule-render.vercel.app/api?type=egg&color=99CCFF&height=50&section=footer"/>
