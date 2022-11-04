import ele
import foyer


def espaloma_to_foyer_xml():
    pass


def mbuild_to_foyer_xml(compound, type_params, bond_params, angle_params, dihedral_params):
    pass


def parmed_to_foyer_xml(structure, ff, file_name):
    """
    """
    # Get needed information from the Parmed structure:
    atom_types = tuple(set(a.type for a in structure.atoms))

    bond_types = set()
    for i in structure.bonds:
        bond_types.add((i.atom1.type, i.atom2.type))
    bond_types = tuple(bond_types)

    angle_types = set()
    for i in structure.angles:
        angle_types.add((i.atom1.type, i.atom2.type, i.atom3.type))
    angle_types = tuple(angle_types)

    dihedral_types = set()
    if "OPLS" in ff.name:
        for i in structure.rb_torsions:
            dihedral_types.add(
                (i.atom1.type, i.atom2.type, i.atom3.type, i.atom4.type)
            )
    dihedral_types = tuple(dihedral_types)
    
    # Write out a new XML file
    with open(file_name, "w") as f:
        f.write(f'<ForceField name="{ff.name}" version="{ff.version}" combining_rule="{ff.combining_rule}">\n')
        f.write("\t<AtomTypes>\n")
        for atom in atom_types:
            atom_type = ff.atomTypeClasses[atom]
            element=ff.atomTypeElements[atom]
            mass = ele.element_from_symbol(element).mass
            _def = ff.atomTypeDefinitions[atom]
            desc = ff.atomTypeDesc[atom]
            line = write_atom_type(
                name=atom,
                atom_type=ff.atomTypeClasses[atom],
                element=ff.atomTypeElements[atom],
                mass=mass,
                _def=_def,
                desc=desc
            )
            f.write(line)
        f.write("\t</AtomTypes>\n")
        # Write out bond parameters
        f.write("\t<HarmonicBondForce>\n")
        for bond in bond_types:
            params = ff.get_parameters("harmonic_bonds", bond)
            class1 = ff.atomTypeClasses[bond[0]]
            class2 = ff.atomTypeClasses[bond[1]]
            line = write_harmonic_bond(class1=class1, class2=class2, l0=params["length"], k=params["k"])
            f.write(line)
        f.write("\t</HarmonicBondForce>\n")
        # Write out angle parameters
        f.write("\t<HarmonicAngleForce>\n")
        for angle in angle_types:
            params=ff.get_parameters("harmonic_angles", angle)
            class1 = ff.atomTypeClasses[angle[0]]
            class2 = ff.atomTypeClasses[angle[1]]
            class3 = ff.atomTypeClasses[angle[2]]
            line = write_harmonic_angle(
                class1=class1, class2=class2, class3=class3, t0=params["theta"], k=params["k"]
            )
            f.write(line)
        f.write("\t</HarmonicAngleForce>\n")
        # Write out dihedral/torsion parameters
        if "OPLS" in ff.name:
            f.write("\t<RBTorsionForce>\n")
            for dihedral in dihedral_types:
                params = ff.get_parameters("rb_propers", dihedral)
                class1 = ff.atomTypeClasses[dihedral[0]]
                class2 = ff.atomTypeClasses[dihedral[1]]
                class3 = ff.atomTypeClasses[dihedral[2]]
                class4 = ff.atomTypeClasses[dihedral[3]]
                line = write_rb_torsion(
                    class1=class1,
                    class2=class2,
                    class3=class3,
                    class4=class4,
                    c0=params["c0"],
                    c1=params["c1"],
                    c2=params["c2"],
                    c3=params["c3"],
                    c4=params["c4"],
                    c5=params["c5"],
                )
                f.write(line)
            f.write("\t</RBTorsionForce>\n")
        # Write out non-bonded parameters 
        f.write(f'\t<NonbondedForce coulomb14scale="{ff.coulomb14scale}" lj14scale="{ff.lj14scale}">\n')
        for atom in atom_types:
            params = ff.get_parameters("atoms", atom)
            line = write_non_bonded(
                name=atom, charge=params["charge"], sigma=params["sigma"], epsilon=params["epsilon"]
            )
            f.write(line)
        f.write('\t</NonbondedForce>\n')
        f.write('</ForceField>')


def write_atom_type(name, atom_type, element, mass, _def=None, desc=None):
    line = f'\t\t<Type name="{name}" class="{atom_type}" element="{element}" mass="{mass}" def="{_def}" desc="{desc}"/>\n'
    return line


def write_harmonic_bond(class1, class2, l0, k):
    line = f'\t\t<Bond class1="{class1}" class2="{class2}" length="{l0}" k="{k}"/>\n'
    return line


def write_harmonic_angle(class1, class2, class3, t0, k):
    line = f'\t\t<Angle class1="{class1}" class2="{class2}" class3="{class3}" angle="{t0}" k="{k}"/>\n'
    return line


def write_non_bonded(name, charge, sigma, epsilon):
    line = f'\t\t<Atom type="{name}" charge="{charge}" sigma="{sigma}" epsilon="{epsilon}"/>\n'
    return line


def write_rb_torsion(class1, class2, class3, class4, c0, c1, c2, c3, c4, c5):
    line = f'\t\t<Proper class1="{class1}" class2="{class2}" class3="{class3}" class4="{class4}" c0="{c0}" c1="{c1}" c2="{c2}" c3="{c3}" c4="{c4}" c5="{c5}"/>\n'
    return line


def write_dihedral():
    pass
