import ele
import foyer
from warnings import warn

def espaloma_to_foyer_xml():
    pass


def parmed_to_foyer_xml(structure, ff, file_name, torsion_type=None):
    """Given a typed Parmed structure, and the Foyer forcefield applied,
    creates and saves a trucated Foyer xml file containing only the
    parameters used in the system.

    Parameters:
    -----------
    structure : pmd.Structure; required
        Parmed structure that has the types, bonds, angles, etc... populated
    ff : Foyer.forcefield.Forcefield; required
        A Foyer Forcefield object that was used when creating the typed Parmed structure
    file_name : str; required
        The file path and name to save the truncated xml file to

    Example:
    --------
    import mbuild as mb
    import foyer
    alkane = mb.load("CCCCC", smiles=True)
    opls = foyer.Forcefield(name="oplsaa")
    alkane_typed = opls.apply(alkane)
    parmed_to_foyer_xml(alkane_typed, opls, "alkane_opls.xml")

    """
    if ff.name and "OPLS" in ff.name:
        if torsion_type and torsion_type != "rb":
            warn(
                    "The forcefield provied appears to be an OPLS type "
                    "which uses the Ryckaert-Bellemans form for dihedrals. "
                    "Setting the torsion type to RB."
            )
        torsion_type = "rb"
    elif ff.name and "GAFF" in ff.name:
        if torsion_type and torsion_type != "periodic":
            warn(
                    "The forcefield provied appears to be a GAFF type "
                    "which uses the periodic form for dihedrals. "
                    "Setting the torsion type to periodic."
            )
        torsion_type = "periodic"

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
    if torsion_type.lower() == "rb":
        for i in structure.rb_torsions:
            dihedral_types.add(
                (i.atom1.type, i.atom2.type, i.atom3.type, i.atom4.type)
            )
    else:
        for i in structure.dihedrals:
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
            try:
                _def = ff.atomTypeDefinitions[atom]
            except KeyError:
                _def = ""
            try:
                desc = ff.atomTypeDesc[atom]
            except KeyError:
                desc = ""
            try:
                override_types = list(ff.atomTypeOverrides[atom])
                override_keep = [i for i in override_types if i in atom_types]
                override = ",".join(override_keep)
            except KeyError:
                override = ""
            line = write_atom_type(
                name=atom,
                atom_type=ff.atomTypeClasses[atom],
                element=ff.atomTypeElements[atom],
                mass=mass,
                _def=_def,
                desc=desc,
                override=override
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
        if torsion_type == "rb":
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
        elif torsion_type == "periodic":
            f.write("\t<PeriodicTorsionForce>\n")
            for dihedral in dihedral_types:
                params = ff.get_parameters("periodic_propers", dihedral)
                if not params:
                    continue
                class1 = ff.atomTypeClasses[dihedral[0]]
                class2 = ff.atomTypeClasses[dihedral[1]]
                class3 = ff.atomTypeClasses[dihedral[2]]
                class4 = ff.atomTypeClasses[dihedral[3]]
                periodicity = params["periodicity"]
                k = params["k"]
                line = write_periodic_dihedral(
                    class1=class1,
                    class2=class2,
                    class3=class3,
                    class4=class4,
                    periodicity=params["periodicity"],
                    k=params["k"],
                    phase=params["phase"]
                )
                f.write(line)
            f.write("\t</PeriodicTorsionForce>\n")

    

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


def mbuild_to_foyer_xml(
    file_name=None,
    compound=None,
    bond_params=None,
    angle_params=None,
    dihedral_params=None,
    dihedral_type="periodic",
    non_bonded_params=None,
    combining_rule="geometric",
    name="",
    version="",
    coulomb14scale=1.0,
    lj14scale=1.0
):
    """
    """
    
    particle_types = tuple(set(p.name for p in compound.particles()))
    particle_masses = []
    for _type in particle_types:
        mass = [p.mass for p in compound.particles_by_name(_type)][0]
        particle_masses.append(mass)
    
        
    with open(file_name, "w") as f:
        f.write(f'<ForceField name="{name}" version="{version}" combining_rule="{combining_rule}">\n')
        f.write("\t<AtomTypes>\n")
        # Write the particle types
        for idx, p in enumerate(particle_types):
            line = write_atom_type(
                name=p,
                atom_type=p,
                element=f"_{p}",
                mass=particle_masses[idx],
                _def=f"_{p}",
            )
            f.write(line)
        f.write("\t</AtomTypes>\n")

        # Write out harmonic bond parameters 
        f.write("<HarmonicBondForce>\n")
        for b in bond_params:
            line = write_harmonic_bond(
                class1=b[0],
                class2=b[1],
                l0=bond_params[b]["l0"],
                k=bond_params[b]["k"]
            )
            f.write(line)
        f.write("</HarmonicBondForce>\n")

        # Write out harmonic angle parameters 
        f.write("<HarmonicAngleForce>\n")
        for a in angle_params:
            line = write_harmonic_angle(
                class1=a[0],
                class2=a[1],
                class3=a[2],
                t0=angle_params[a]["t0"],
                k=angle_params[a]["k"]
            )
            f.write(line)
        f.write("</HarmonicAngleForce>\n")

        # Write out dihedral parameters 
        if dihedral_type == "periodic":
            f.write("<PeriodicTorsionForce>\n")
            for d in dihedral_params:
                line = write_periodic_dihedral(
                    class1=d[0],
                    class2=d[1],
                    class3=d[2],
                    class4=d[3],
                    periodicity=dihedral_params[d]["periodicity"],
                    k=dihedral_params[d]["k"],
                    phase=dihedral_params[d]["phase"],
                )
                f.write(line)
            f.write("</PeriodicTorsionForce>\n")

        # Write out non-bonded parameters 
        f.write(f'\t<NonbondedForce coulomb14scale="{coulomb14scale}" lj14scale="{lj14scale}">\n')
        for a in non_bonded_params:
            line = write_non_bonded(
                name=a,
                charge=non_bonded_params[a]["charge"],
                sigma=non_bonded_params[a]["sigma"],
                epsilon=non_bonded_params[a]["epsilon"]
            )
            f.write(line)
        f.write('\t</NonbondedForce>\n')
        f.write('</ForceField>')


def write_atom_type(name, atom_type, element, mass, _def="", desc="", override=""):
    """Creates a line for an atom type following the foyer-xml format"""
    line = f'\t\t<Type name="{name}" class="{atom_type}" element="{element}" mass="{mass}" def="{_def}" desc="{desc}" overrides="{override}"/>\n'
    return line


def write_harmonic_bond(class1, class2, l0, k):
    """
    """
    line = f'\t\t<Bond class1="{class1}" class2="{class2}" length="{l0}" k="{k}"/>\n'
    return line


def write_harmonic_angle(class1, class2, class3, t0, k):
    """
    """
    line = f'\t\t<Angle class1="{class1}" class2="{class2}" class3="{class3}" angle="{t0}" k="{k}"/>\n'
    return line


def write_non_bonded(name, charge, sigma, epsilon):
    """
    """
    line = f'\t\t<Atom type="{name}" charge="{charge}" sigma="{sigma}" epsilon="{epsilon}"/>\n'
    return line


def write_rb_torsion(class1, class2, class3, class4, c0, c1, c2, c3, c4, c5):
    """
    """
    line = f'\t\t<Proper class1="{class1}" class2="{class2}" class3="{class3}" class4="{class4}" c0="{c0}" c1="{c1}" c2="{c2}" c3="{c3}" c4="{c4}" c5="{c5}"/>\n'
    return line


def write_periodic_dihedral(class1, class2, class3, class4, periodicity, k, phase):
    """
    """
    # Extend lists of periodicity, k and phase to ensure they are always len 4
    periodicity.extend([0] * (4-len(periodicity)))
    k.extend([0] * (4-len(k)))
    phase.extend([0] * (4-len(phase)))
    line = f'\t\t<Proper class1="{class1}" class2="{class2}" class3="{class3}" class4="{class4}" periodicity1="{periodicity[0]}" k1="{k[0]}" phase1="{phase[0]}" periodicity2="{periodicity[1]}" k2="{k[1]}" phase2="{phase[1]}" periodicity3="{periodicity[2]}" k3="{k[2]}" phase3="{phase[2]}" periodicity4="{periodicity[3]}" k4="{k[3]}" phase4="{phase[3]}"/>\n'
    return line
