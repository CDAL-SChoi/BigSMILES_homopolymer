from setuptools import setup, find_packages

setup(
    name='BigSMILES_homopolymer',
    version='1.0.0',
    description='An automated BigSMILES conversion workflow for homopolymeric macromolecules',
    url='https://github.com/CDAL-SChoi/BigSMILES_homopolymer',
    author='Sunho Choi',
    author_email='schoi_@korea.ac.kr',
    license='schoi',
    packages=find_packages(),
    install_requires=['numpy', 'openpyxl==3.0.9', 'pandas', 'rdkit==2022.9.5', 'regex==2022.6.2'],
)
