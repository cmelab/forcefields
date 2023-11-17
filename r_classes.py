import mbuild as mb
import warnings 
import numpy
import ipywidgets as widgets
warnings.filterwarnings('ignore')

class c11_bo(mb.Compound):
    def __init__(self):
        super(c11_bo, self).__init__()
        self.add(mb.load("CCCCCCCCCCCC(CCCC)CCCCCC",smiles=True))
        carbon1 = list(self.particles_by_name('C'))[0]
        self.add(mb.Port(anchor=carbon1, orientation=[1, -1, -1], separation=0.08), "p1")
        H1 = list(self.particles_by_name('H'))[0]
        self.remove(H1)

class HD(mb.Compound):
    def __init__(self):
        super(HD, self).__init__()
        self.add(mb.load("CC(CCCCCC)CCCCCCCC",smiles=True))
        carbon1 = list(self.particles_by_name('C'))[0]
        self.add(mb.Port(anchor=carbon1, orientation=[1, 1, 1], separation=0.08), "p1")
        H1 = list(self.particles_by_name('H'))[0]
        self.remove(H1)

class ODD(mb.Compound):
    def __init__(self):
        super(ODD, self).__init__()
        self.add(mb.load("CC(CCCCCCCCCC)CCCCCCCC",smiles=True))
        carbon1 = list(self.particles_by_name('C'))[0]
        self.add(mb.Port(anchor=carbon1, orientation=[-1, -1, -1], separation=0.08), "p1")
        H1 = list(self.particles_by_name('H'))[0]
        self.remove(H1)

class C1BO(mb.Compound):
    def __init__(self):
        super(C1BO, self).__init__()
        self.add(mb.load("CC(CCCC)CCCCCC",smiles=True))
        carbon1 = list(self.particles_by_name('C'))[0]
        self.add(mb.Port(anchor=carbon1, orientation=[-1, -1, -1], separation=0.08), "p1")
        H1 = list(self.particles_by_name('H'))[0]
        self.remove(H1)

class C3BO(mb.Compound):
    def __init__(self):
        super(C3BO, self).__init__()
        self.add(mb.load("CCCC(CCCC)CCCCCC",smiles=True))
        carbon1 = list(self.particles_by_name('C'))[0]
        self.add(mb.Port(anchor=carbon1, orientation=[-1, -1, -1], separation=0.08), "p1")
        H1 = list(self.particles_by_name('H'))[0]
        self.remove(H1)


class C4BO(mb.Compound):
    def __init__(self):
        super(C4BO, self).__init__()
        self.add(mb.load("CCCCC(CCCC)CCCCCC",smiles=True))
        carbon1 = list(self.particles_by_name('C'))[0]
        self.add(mb.Port(anchor=carbon1, orientation=[-1, -1, -1], separation=0.08), "p1")
        H1 = list(self.particles_by_name('H'))[0]
        self.remove(H1)

class C5BO(mb.Compound):
    def __init__(self):
        super(C5BO, self).__init__()
        self.add(mb.load("CCCCCC(CCCC)CCCCCC",smiles=True))
        carbon1 = list(self.particles_by_name('C'))[0]
        self.add(mb.Port(anchor=carbon1, orientation=[1, 1, 1], separation=0.08), "p1")
        H1 = list(self.particles_by_name('H'))[0]
        self.remove(H1)



class C16(mb.Compound):
    def __init__(self):
        super(C16, self).__init__()
        self.add(mb.load("CCCCCCCCCCCCCCCC",smiles=True))
        carbon1 = list(self.particles_by_name('C'))[0]
        self.add(mb.Port(anchor=carbon1, orientation=[1, 1, 1], separation=0.08), "p1")
        H1 = list(self.particles_by_name('H'))[0]
        self.remove(H1)


class ene_ODD(mb.Compound):
    def __init__(self):
        super(ene_ODD, self).__init__()
        self.add(mb.load("CCCCCCCCCCC(CCCCCCCC)CCCC(CCCCCCCC)CCCCCCCCCC",smiles=True))
        carbon1 = self[20]
        self.add(mb.Port(anchor=carbon1, orientation=[1, -1, 1], separation=0.08), "p1")
        self.add(mb.Port(anchor=carbon1, orientation=[1, -2, 1], separation=0.08), "p2")
        H1 = self[82]
        H2 = self[83]
        self.remove(H1)
        self.remove(H2)
        
        
class ene_HD(mb.Compound):
    def __init__(self):
        super(ene_HD, self).__init__()
        self.add(mb.load("CCCCCCCCC(CCCCCC)CCCC(CCCCCC)CCCCCCCC",smiles=True))
        carbon1 = self[16]
        self.add(mb.Port(anchor=carbon1, orientation=[0, 1, -1], separation=0.08), "p1")
        self.add(mb.Port(anchor=carbon1, orientation=[-1, -1, -1], separation=0.08), "p2")
        H1 = self[66]
        H2 = self[67]
        self.remove(H1)
        self.remove(H2)
