import warnings 
warnings.filterwarnings('ignore')
import mbuild as mb
import numpy
from mbuild.lib.recipes.polymer import Polymer
import espaloma as esp
from writers import foyer_xml_writer
from writers.foyer_xml_writer import mbuild_to_foyer_xml
from functions import bond_walk
from functions.bond_walk import MadAtom, MadBond, BondWalker
import ipywidgets as widgets
import os
import torch
import parmed as pmd
import networkx  as nx
if not os.path.exists("espaloma_model.pt"):
    os.system("wget http://data.wangyq.net/espaloma_model.pt")

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
        chain.energy_minimize()
    return chain

#def get_aromatic_ids(compound):
#    aromatic_particle_indices = set()
#    for idx, p in enumerate(compound.particles()):
#        for bond in compound.bonds(return_bond_order=True):
#            if p in [bond[0], bond[1]] and bond[2]["bond_order"] == "aromatic":
#                aromatic_particle_indices.add(idx)
#    return aromatic_particle_indices
#
#
#def fix_rdkit_aromaticity(rdmol, particle_indices):
#    for idx, a in enumerate(rdmol.GetAtoms()):
#        if idx in particle_indices:
#            a.SetIsAromatic(True)

def espaloma(MONOMER,XML_FILEPATH,TYPED_FILEPATH,DIMER):
    '''if you plan to parameterize a polymer made up of your monomer set dimer = True, 
    if you are simulating only the monomer set dimer = False
    '''
    if DIMER == True:
        #mb_mol = build_chain(monomer=MONOMER,length=2,min_energy=False)
        monomer = build_chain(monomer=MONOMER,length=1,min_energy=False)
        dimer = build_chain(monomer=MONOMER,length=2,min_energy=False)
        mb_mol = mb.Compound()
        mb_mol.add([monomer,dimer])
        monomer.translate([1,1,1])
        dimer.translate([-1,-1,-1])
    else:
        mb_mol = MONOMER
    rdkit_mol = mb_mol.to_rdkit()
    #aromatic_indices =  get_aromatic_ids(mb_mol)
    #fix_rdkit_aromaticity(rdkit_mol,aromatic_indices)
    from openff.toolkit.topology import Molecule
    comp = Molecule.from_rdkit(rdkit_mol,allow_undefined_stereo=True,hydrogens_are_explicit=True)
    bonds = [b for b in comp.bonds]
    for i in range(len(bonds)):
        bonds[i].bond_order = 1
    b= BondWalker(comp)
    molecule = b.fill_in_bonds()
    molecule_graph = esp.Graph(molecule)
    espaloma_model = esp.get_model("latest")
    espaloma_model(molecule_graph.heterograph)
    openmm_system = esp.graphs.deploy.openmm_system_from_graph(molecule_graph,charge_method="nn")
    
    # Store the results for each in something more accessible
    pair_forces = openmm_system.getForces()[0]
    angle_forces = openmm_system.getForces()[2]
    bond_forces = openmm_system.getForces()[3]
    torsion_forces = openmm_system.getForces()[1]
    
    # get a parmed structure from openmm 
    topology = molecule.to_topology()
    openmm_topology = topology.to_openmm()
    
    structure = pmd.openmm.load_topology(topology=openmm_topology, system=openmm_system)
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
            
    #generate type_map:
    Gopenmm = nx.Graph()
    Gparmed = nx.Graph()
    
    #openmm:
    for i in range(bond_forces.getNumBonds()):
        Gopenmm.add_edge(bond_forces.getBondParameters(index=i)[0],bond_forces.getBondParameters(index=i)[1])
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
        for i in range(pair_forces.getNumParticles()):
            pair_parms = pair_forces.getParticleParameters(index=i)
            sigma = pair_parms[1]/pair_parms[1].unit
            epsilon = pair_parms[2]/pair_parms[2].unit
            if (sigma, epsilon) not in particle_types: 
                particle_types.append((sigma, epsilon))
            type_map[molecule.atoms[i].molecule_atom_index] = "".join([molecule.atoms[i].name , 
                                                                       str(particle_types.index((sigma, epsilon)))])
    
    
    #generate dictionaries:
    
    bond_types = []
    bond_dict = dict() 
    
    for i in range(bond_forces.getNumBonds()):
        bond_parms = bond_forces.getBondParameters(index=i)
        l0 = bond_parms[2]/bond_parms[2].unit
        k = bond_parms[3]/bond_parms[3].unit
        bond_dict[type_map[bond_parms[0]],type_map[bond_parms[1]]] = {'k':k,'l0':l0}
    
    
    angle_types = []
    angle_dict = dict()
    
    for i in range(angle_forces.getNumAngles()):
        angle_parms = angle_forces.getAngleParameters(index=i)
        k = angle_parms[4]/angle_parms[4].unit
        t0 = angle_parms[3]/angle_parms[3].unit  
        angle_dict[type_map[angle_parms[0]],type_map[angle_parms[1]],type_map[angle_parms[2]]] = {'k':k,'t0':t0}
    
    
    
    dihedral_types = []
    dihedral_dict = {}
    
    for i in range(torsion_forces.getNumTorsions()):
        if i%6==0:
            periodicity=[]
            phase = []
            k = []
        dihedral_parms = torsion_forces.getTorsionParameters(index=i)
        periodicity.append(dihedral_parms[4])  
        phase.append( dihedral_parms[5]/dihedral_parms[5].unit)
        k.append(dihedral_parms[6]/dihedral_parms[6].unit)
        dt = (type_map[dihedral_parms[0]],type_map[dihedral_parms[1]],type_map[dihedral_parms[2]],
                      type_map[dihedral_parms[3]])
       
    
        if periodicity[-1]==6:
            dihedral_dict[dt] = {'periodicity':periodicity,'k':k,'phase':phase}
    
    
    nonbonded_types = []
    nonbonded_dict = {}
    
    for i in range(pair_forces.getNumParticles()):
        nonbonded_parms = pair_forces.getParticleParameters(index=i)
        charge = nonbonded_parms[0]/nonbonded_parms[0].unit
        sigma = nonbonded_parms[1]/nonbonded_parms[1].unit
        epsilon = nonbonded_parms[2]/nonbonded_parms[2].unit
        nonbonded_types.append((charge,sigma,epsilon))
        nonbonded_dict[(type_map[i])]={'charge':charge,'sigma':sigma,'epsilon':epsilon}
    
    for index in type_map:
        #print(index, type_map[index],mb_mol[index].name)
        mb_mol[index].name = type_map[index]


    if DIMER == True and mb_mol.children[1].n_particles > mb_mol.children[0].n_particles:
        dim = mb_mol.children[1]
        mb_mol.remove(objs_to_remove=dim)
    
    mbuild_to_foyer_xml(
        file_name=XML_FILEPATH,
        compound=mb_mol,
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
    
    mb_mol.save(TYPED_FILEPATH,overwrite=True)

    if os.path.exists(TYPED_FILEPATH):
        print('Typed mol2 file created.')
    else: 
        print('Failed to create typed mol2')
    if os.path.exists(XML_FILEPATH):
        print('xml file created.')
    else: 
        print('Failed to create xml file')



def build_polymer(monomer, length,bond_indices, separation,replace,orientations,min_energy):
    chain = Polymer()
    chain.add_monomer(compound=monomer,
                 indices=bond_indices,
                 separation=separation,
                 replace=replace,
                 orientation=orientations)
    chain.build(n=length)
    if chain[-2].name == 'H' and chain[-1].name == 'H':
        chain[-2].name = monomer[bond_indices[0]].name
        chain[-1].name = monomer[bond_indices[1]].name
    if min_energy == True:
        chain.energy_minimize()
    return chain
