<img src="https://capsule-render.vercel.app/api?type=egg&color=99CCFF&height=50&section=header"/>


# An automated BigSMILES conversion workflow and dataset for homopolymeric macromolecules

This repository contains the code of generating BigSMILES from the SMILES representations of homopolymers, and the opposite process as well.
The dataset is available from the [Figshare](https://figshare.com/s/813ca7794bd9257e9843).

**Abstract:**
SMILES has been utilized in various artificial intelligence analyses due to its ability to represent chemical structures in line notation. However, its ease of representation is limited, particularly for macromolecules, which has led to the proposal of BigSMILES as a representation method suitable for macromolecules. Nevertheless, research in BigSMILES has been limited so far due to the requirement for manual conversion. In this paper, we present a conversion workflow of BigSMILES, focusing on its automated generation from SMILES representations of homopolymers, and provide the BigSMILES representations for all 4,927,181 records enabling immediate use for a variety of research and development application. Our study introduces a meticulous validation process to ensure the accuracy, interchangeability, and robustness of the conversion. Additionally, we provide a systematic overview of the code and functions utilized, emphasizing their relevance in the context of BigSMILES generation. We anticipate that this advancement would significantly aid researchers, facilitating further studies in BigSMILES representation, including potential applications in deep learning and further extension to complex structures such as copolymers.

## Install BigSMILES_homopolymer package
'''py
pip install git+https://github.com/CDAL-SChoi/BigSMILES_homopolymer.git
'''


## The workflow of the s2bigs.py code
![figure 3](https://github.com/CDAL-SChoi/BigSMILES_homopolymer/assets/50295574/ec6b6b03-1387-459f-b823-ba429bd91ebd)



<img src="https://capsule-render.vercel.app/api?type=egg&color=99CCFF&height=50&section=footer"/>
