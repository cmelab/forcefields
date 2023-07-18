import warnings 
warnings.filterwarnings('ignore')
import mbuild as mb
import numpy
from frag_classes import IDT, CPDT, DPP, BT, PT, FBT, thiophene, pyridine, CPDT_eneHD 
from r_classes import c11_bo,HD,ODD,C1BO,C3BO,C4BO,C5BO, C16
from mon_classes import PCPDT_PT_eneHD,PCPDTPT_eneODD,PCPDTPT_nC16,PIDTBT_nC16, PCPDTPT_ODD,PIDTCPDT_C11BO, P3HT, PCPDTPT_HD, PDPPPyT_ODD, PCPDTFBT_C1_BO, PCPDTFBT_C3_BO, PCPDTFBT_C4_BO, PCPDTFBT_C5_BO, PIDTFBT_C11_BO, PCPDTFBT_C11_BO, BDT_TPD, perylene
from mbuild.lib.recipes.polymer import Polymer

from writers import foyer_xml_writer
from writers.foyer_xml_writer import mbuild_to_foyer_xml
import ipywidgets as widgets
import os
import torch
import parmed as pmd
import networkx  as nx
if not os.path.exists("espaloma_model.pt"):
    os.system("wget http://data.wangyq.net/espaloma_model.pt")
    
    
from bondwalk import bond_walk
from bondwalk.bond_walk import MadAtom, MadBond, BondWalker
import espaloma as esp
from openff.toolkit.topology import Molecule


def build_chain(fragment,length,min_energy):
        chain = Polymer()
        chain.add_monomer(compound=fragment,
                     indices=fragment.bond_indices,
                     separation=fragment.separation,
                     replace=fragment.replace,
                     orientation=fragment.orientations)
        chain.build(n=length)
        if min_energy == True:
            chain.energy_minimize()
        return chain



polymer = build_chain(fragment=P3HT(),
	length = 2,
	min_energy = True)
