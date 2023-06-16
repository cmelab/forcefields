from writers import foyer_xml_writer
from writers.foyer_xml_writer import parmed_to_foyer_xml, mbuild_to_foyer_xml
from bondwalk import bond_walk
from bondwalk.bond_walk import MadAtom, MadBond, BondWalker

import warnings
warnings.filterwarnings("ignore")
import ele
import espaloma as esp
import forcefield_utilities as ffutils
import foyer
import gmso
import mbuild as mb
from mbuild.lib.recipes import Polymer
from mbuild.formats.hoomd_forcefield import create_hoomd_forcefield
import numpy as np
import torch
from openff.toolkit.topology import Molecule
from mbuild.formats.hoomd_forcefield import create_hoomd_forcefield
import hoomd
import gsd.hoomd
import matplotlib.pyplot as plt
import rdkit
from rdkit import Chem

import os

if not os.path.exists("espaloma_model.pt"):
    os.system("wget http://data.wangyq.net/espaloma_model.pt")


def openff_Molecule(MOL2FILEPATH,PDBFILEPATH):
    rdmol = Chem.MolFromMol2File(MOL2FILEPATH)
    from_rdmol = Molecule.from_rdkit(rdmol)
    atom_list = [a for a in from_rdmol.atoms if 1 not in [a.atomic_number]]
    for atom in atom_list:
        atom.formal_charge = 0
    b = BondWalker(from_rdmol)
    comp = b.fill_in_bonds()
    from_rdmol.to_file(PDBFILEPATH,file_format="pdb")
    return comp

def espaloma(molecule):
    molecule_graph = esp.Graph(molecule)
    espaloma_model = torch.load("espaloma_model.pt")
    espaloma_model(molecule_graph.heterograph)
    openmm_system = esp.graphs.deploy.openmm_system_from_graph(molecule_graph,charge_method="nn")
    return openmm_system


def dictionaries(typemap,BondForces,AngleForces,TorsionForces,PairForces):
    bond_types = []
    bond_dict = dict() 
    for i in range(BondForces.getNumBonds()):
        bond_parms = BondForces.getBondParameters(index=i)
        l0 = bond_parms[2]/bond_parms[2].unit
        k = bond_parms[3]/bond_parms[3].unit
        bond_dict[typemap[bond_parms[0]],typemap[bond_parms[1]]] = {'k':k,'l0':l0}
        
    angle_types = []
    angle_dict = dict()
    for i in range(AngleForces.getNumAngles()):
        angle_parms = AngleForces.getAngleParameters(index=i)
        k = angle_parms[4]/angle_parms[4].unit
        t0 = angle_parms[3]/angle_parms[3].unit  
        angle_dict[typemap[angle_parms[0]],typemap[angle_parms[1]],typemap[angle_parms[2]]] = {'k':k,'t0':t0}
    
    dihedral_types = []
    dihedral_dict = {}
    
    for i in range(TorsionForces.getNumTorsions()):
        if i%6==0:
            periodicity=[]
            phase = []
            k = []
        dihedral_parms = TorsionForces.getTorsionParameters(index=i)
        periodicity.append(dihedral_parms[4])  
        phase.append( dihedral_parms[5]/dihedral_parms[5].unit)
        k.append(dihedral_parms[6]/dihedral_parms[6].unit)
        dt = (typemap[dihedral_parms[0]],typemap[dihedral_parms[1]],typemap[dihedral_parms[2]],
                      typemap[dihedral_parms[3]])
       
    
        if periodicity[-1]==6:
            dihedral_dict[dt] = {'periodicity':periodicity,'k':k,'phase':phase}
    
    
    nonbonded_types = []
    nonbonded_dict = {}
    
    for i in range(PairForces.getNumParticles()):
        nonbonded_parms = PairForces.getParticleParameters(index=i)
        charge = nonbonded_parms[0]/nonbonded_parms[0].unit
        sigma = nonbonded_parms[1]/nonbonded_parms[1].unit
        epsilon = nonbonded_parms[2]/nonbonded_parms[2].unit
        nonbonded_types.append((charge,sigma,epsilon))
        nonbonded_dict[(typemap[i])]={'charge':charge,'sigma':sigma,'epsilon':epsilon}

    return bond_dict,angle_dict,dihedral_dict,nonbonded_dict

def typing(BondForces,PairForces,molecule,omm_system):
    import parmed as pmd
    topology = molecule.to_topology()
    openmm_topology = topology.to_openmm()
    structure = pmd.openmm.load_topology(topology=openmm_topology, system=omm_system)
    structure.bonds.sort(key=lambda x: x.atom1.idx)
    
    
    for i in range(len(molecule.atoms)):
        if molecule.atoms[i].atomic_number == 6:
            molecule.atoms[i].name = 'C'
        if molecule.atoms[i].atomic_number == 1:
            molecule.atoms[i].name = 'H'
        if molecule.atoms[i].atomic_number == 7:
            molecule.atoms[i].name = 'N'
        if molecule.atoms[i].atomic_number == 16:
            molecule.atoms[i].name = 'S'
        if molecule.atoms[i].atomic_number == 8:
            molecule.atoms[i].name = 'O'
        if molecule.atoms[i].atomic_number == 9:
            molecule.atoms[i].name = 'F'
    
    import networkx  as nx
    Gopenmm = nx.Graph()
    Gparmed = nx.Graph()
    #openmm:
    for i in range(BondForces.getNumBonds()):
        Gopenmm.add_edge(BondForces.getBondParameters(index=i)[0],BondForces.getBondParameters(index=i)[1])
    #parmed
    for b in structure.bonds:
        Gparmed.add_edge(b.atom1.idx,b.atom2.idx)
        
    particle_types = []
    type_map = dict()
    
    #nx.rooted_tree_isomorphism
    #in here we still need to check that one known index on one corresponds to the same index on the other....
    tree_openmm = nx.bfs_tree(Gopenmm,0)
    tree_parmed = nx.bfs_tree(Gparmed,0)
    if nx.is_isomorphic(Gopenmm,Gparmed):
    #if nx.isomorphism.tree_isomorphism(tree_openmm,tree_parmed):  <- want this work
        for i in range(PairForces.getNumParticles()):
            pair_parms = PairForces.getParticleParameters(index=i)
            sigma = pair_parms[1]/pair_parms[1].unit
            epsilon = pair_parms[2]/pair_parms[2].unit
            if (sigma, epsilon) not in particle_types: 
                particle_types.append((sigma, epsilon))
            type_map[molecule.atoms[i].molecule_atom_index] = "".join([molecule.atoms[i].name , 
                                                                       str(particle_types.index((sigma, epsilon)))])
    return type_map

def esp_to_xml(MOL2FILEPATH,XMLFILEPATH,PDBFILEPATH,TYPEDFILEPATH):
    off_molecule = openff_Molecule(MOL2FILEPATH=MOL2FILEPATH,PDBFILEPATH=PDBFILEPATH)
    openmm_system = espaloma(molecule=off_molecule)
    pair_forces = openmm_system.getForces()[1]
    angle_forces = openmm_system.getForces()[3]
    bond_forces = openmm_system.getForces()[2]
    torsion_forces = openmm_system.getForces()[0]
    type_map = typing(BondForces=bond_forces,PairForces=pair_forces,molecule=off_molecule,omm_system=openmm_system)
    comp_rename = mb.load(PDBFILEPATH)
    for index in type_map:
        comp_rename[index].name = type_map[index]
    bond_dict = dictionaries(typemap=type_map,BondForces=bond_forces,AngleForces=angle_forces,
                             TorsionForces=torsion_forces,PairForces=pair_forces)[0]
    angle_dict = dictionaries(typemap=type_map,BondForces=bond_forces,AngleForces=angle_forces,
                              TorsionForces=torsion_forces,PairForces=pair_forces)[1]
    dihedral_dict = dictionaries(typemap=type_map,BondForces=bond_forces,AngleForces=angle_forces,
                                 TorsionForces=torsion_forces,PairForces=pair_forces)[2]
    nonbonded_dict = dictionaries(typemap=type_map,BondForces=bond_forces,AngleForces=angle_forces,
                                  TorsionForces=torsion_forces,PairForces=pair_forces)[3]
    
    mbuild_to_foyer_xml(
        file_name=XMLFILEPATH, #change this to whatever you want to save your xml file as
        compound=comp_rename,
        bond_params=bond_dict,
        angle_params=angle_dict,
        dihedral_params=dihedral_dict,
        dihedral_type="periodic",
        non_bonded_params=nonbonded_dict,
        combining_rule="geometric",
        name="",
        version="",
        coulomb14scale=1.0,
        lj14scale=1.0)
    comp_rename.save(TYPEDFILEPATH, overwrite=True) #change this to match your molecule name. 
