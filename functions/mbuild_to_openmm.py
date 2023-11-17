import mbuild as mb
from openff.toolkit.topology import Molecule
from bondwalk import bond_walk
from bondwalk.bond_walk import MadAtom, MadBond, BondWalker


def mbuild_to_openmm(SMILES,sdf_filepath):
    comp = Molecule.from_smiles(SMILES)
    bonds = [b for b in comp.bonds]
    for i in range(len(bonds)):
        bonds[i].bond_order = 1
    b= BondWalker(comp)
    molecule = b.fill_in_bonds()
    molecule.to_file(sdf_filepath,file_format='SDF')
