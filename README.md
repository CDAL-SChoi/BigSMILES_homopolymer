<img src="https://capsule-render.vercel.app/api?type=egg&color=99CCFF&height=50&section=header"/>


# Automated BigSMILES conversion workflow and dataset for homopolymeric macromolecules

This repository contains the code of generating BigSMILES from the SMILES representations of homopolymers, and the opposite process as well.
The dataset is available from the [Figshare](https://figshare.com/s/813ca7794bd9257e9843).

Note that this SMILES to BigSMILES converter operates to convert to the BigSMILES version 1.0 notation, as defined in original BigSMILES paper(2019).
To check the detailed BigSMILES line notation rules, visit [documentation](https://olsenlabmit.github.io/BigSMILES/docs/line_notation.html#the-bigsmiles-line-notation) from the BigSMILES development team.

**Abstract:**
The simplified molecular-input line-entry system (SMILES) has been utilized in a variety of artificial intelligence analyses owing to its capability of representing chemical structures using line notation. However, its ease of representation is limited, which has led to the proposal of BigSMILES as an alternative method suitable for the representation of macromolecules. Nevertheless, research on BigSMILES remains limited due to its preprocessing requirements. Thus, this study proposes a conversion workflow of BigSMILES, focusing on its automated generation from SMILES representations of homopolymers. BigSMILES representations for 4,927,181 records are provided, thereby enabling its immediate use for various research and development applications. Our study presents detailed descriptions on a validation process to ensure the accuracy, interchangeability, and robustness of the conversion. Additionally, a systematic overview of utilized codes and functions that emphasizes their relevance in the context of BigSMILES generation are produced. This advancement is anticipated to significantly aid researchers and facilitate further studies in BigSMILES representation, including potential applications in deep learning and further extension to complex structures such as copolymers.

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
![workflow_revised](https://github.com/CDAL-SChoi/BigSMILES_homopolymer/assets/50295574/37062d6c-9489-4115-b0aa-2b2c0627f2f3)




<img src="https://capsule-render.vercel.app/api?type=egg&color=99CCFF&height=50&section=footer"/>
