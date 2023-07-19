import warnings 
warnings.filterwarnings('ignore')
import mbuild as mb
import numpy
from frag_classes import IDT, CPDT, DPP, BT, PT, FBT, thiophene, pyridine, CPDT_eneHD 
from r_classes import c11_bo,HD,ODD,C1BO,C3BO,C4BO,C5BO, C16
from mon_classes import PCPDTPT_nC16,PIDTBT_nC16, PCPDTPT_ODD,PIDTCPDT_C11BO, P3HT, PCPDTPT_HD, PDPPPyT_ODD, PCPDTFBT_C1_BO, PCPDTFBT_C3_BO, PCPDTFBT_C4_BO, PCPDTFBT_C5_BO, PIDTFBT_C11_BO,PCPDT_PT_eneHD, PCPDTFBT_C11_BO, BDT_TPD, perylene
from mbuild.lib.recipes.polymer import Polymer

def build_chain(monomer, length, min_energy):
    chain = Polymer()
    chain.add_monomer(compound=monomer,
                 indices=monomer.bond_indices,
                 separation=monomer.separation,
                 replace=monomer.replace,
                 orientation=monomer.orientations)
    chain.build(n=length)
    if min_energy == True:
        chain.energy_minimize()
    return chain

def mon_to_smiles(fragment):
    monomer = build_chain(fragment,1,min_energy=True)
    dimer = build_chain(fragment,2,min_energy=True)
    
    mon_dim = mb.Compound()
    mon_dim.add([monomer,dimer])
    monomer.translate([3,3,3])
    dimer.translate([-3,-3,-3])
    smiles_string = mon_dim.to_smiles()
    return mon_dim, smiles_string
