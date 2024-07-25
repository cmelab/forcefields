# forcefields

This repo is designed to aid in the proccess of creating forcefield files from scratch using the Espaloma Program written by the Chodera Lab (https://github.com/choderalab/espaloma). 

# Installing Environment 

## For Apple Silicon Mac (ARM):

```
CONDA_SUBDIR=osx-64 conda create -n ff_generation python
```
```
conda activate ff_generation
conda env config vars set CONDA_SUBDIR=osx-64
conda deactivate
conda activate ff_generation
```
```
conda install 'espaloma=0.3.2' mbuild "python=3.10.14" foyer openbabel py3dmol
```
```
pip install notebook
```

# Example 

from functions import MonToSmiles
from functions.MonToSmiles import mon_to_smiles

#you can either import a preexisting monomer class:
from monclasses.py import P3HT

#or you can build the monomer class:

class P3HT(mb.Compound):
    def __init__(self):
        super(P3HT,self).__init__()
        self.add(mb.load("CCCCCCC1=C(SC(=C1))",smiles=True))
        self.bond_indices = [24,25]
        self.orientations = [[0,0,1],[0,0,-1]]
        self.separation = 0.14
        self.replace = True

#The monomer class must include the atomic indices where the polymerizing bonds will be formed as well as the orientations of the #bonds, the bond separation, and whether you are replacing the atoms where the bond will be formed or not. 

smiles = mon_to_smiles(fragment=P3HT())[1]   #this line must be run before importing the espaloma function. 

#The mon_to_smiles function returns both the mbuild compound including both the monomer and dimer as well as the smiles string of the monomer and dimer 

from functions import EspalomaFxn
from functions.EspalomaFxn import espaloma

espaloma(SMILES=smiles,  #from mon_to_smiles function
        XML_FILEPATH='INSERT DESIRED XML FILEPATH',
        TYPED_FILEPATH='INSERT DESIRED MOL2 FILEPATH')
