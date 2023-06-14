import mbuild as mb
import warnings 
import numpy
import ipywidgets as widgets
warnings.filterwarnings('ignore')

class CPDT(mb.Compound):
    def __init__(self):
        super(CPDT, self).__init__()
        self.add(mb.load("C1C3=C(SC=C3)C2=C1C=CS2", smiles=True))
        carbon1 = list(self.particles_by_name('C'))[8]
        carbon2 = list(self.particles_by_name('C'))[3]
        carbon3 = list(self.particles_by_name('C'))[0]
        self.add(mb.Port(anchor=carbon1, orientation=[10, 1, 1], separation=0.08), "p1")
        self.add(mb.Port(anchor=carbon2, orientation=[-10, -1, -1], separation=0.08), "p2")
        self.add(mb.Port(anchor=carbon3, orientation=[-1, 10, -10], separation=0.08), "p3")
        self.add(mb.Port(anchor=carbon3, orientation=[-1, -10, -5], separation=0.08), "p4")
        hydrogen1 = list(self.particles_by_name('H'))[5]
        hydrogen2 = list(self.particles_by_name('H'))[2]
        hydrogen3 = list(self.particles_by_name('H'))[0]
        hydrogen4 = list(self.particles_by_name('H'))[1]
        self.remove(hydrogen1)
        self.remove(hydrogen2)
        self.remove(hydrogen3)
        self.remove(hydrogen4)


class PT(mb.Compound):
    def __init__(self):
        super(PT, self).__init__()
        self.add(mb.load("n1ccc2nsnc2c1",smiles=True))
        carbon1 = list(self.particles_by_name('C'))[4]
        carbon2 = list(self.particles_by_name('C'))[1]
        self.add(mb.Port(anchor=carbon1, orientation=[0,1,0], separation=0.08), "p1")
        self.add(mb.Port(anchor=carbon2, orientation=[0,-1,0], separation=0.08), "p2")
        H1 = list(self.particles_by_name('H'))[2]
        H2 = list(self.particles_by_name('H'))[1]
        self.remove(H1)
        self.remove(H2)

class IDT(mb.Compound):
    def __init__(self):
        super(IDT, self).__init__()
        self.add(mb.load("C1C2=C(C3=CC4=C(C=C31)C5=C(C4)C=CS5)SC=C2", smiles=True))
        carbon1 = list(self.particles_by_name('C'))[14]
        carbon2 = list(self.particles_by_name('C'))[13]
        carbon3 = list(self.particles_by_name('C'))[0]
        carbon4 = list(self.particles_by_name('C'))[0]
        carbon5 = list(self.particles_by_name('C'))[11]
        carbon6 = list(self.particles_by_name('C'))[11]
        self.add(mb.Port(anchor=carbon1, orientation=[10, 1, 1], separation=0.08), "p1")
        self.add(mb.Port(anchor=carbon2, orientation=[-10, -1, -1], separation=0.08), "p2")
        self.add(mb.Port(anchor=carbon3, orientation=[0, -1, 10], separation=0.08), "p3")
        self.add(mb.Port(anchor=carbon4, orientation=[0, -11, -10], separation=0.08), "p4")
        self.add(mb.Port(anchor=carbon5, orientation=[0, 1, -10], separation=0.08), "p5")
        self.add(mb.Port(anchor=carbon6, orientation=[0, 1, 1], separation=0.08), "p6")
        H1 = list(self.particles_by_name('H'))[8]
        H2 = list(self.particles_by_name('H'))[7]
        H3 = list(self.particles_by_name('H'))[0]
        H4 = list(self.particles_by_name('H'))[1]
        H5 = list(self.particles_by_name('H'))[4]
        H6 = list(self.particles_by_name('H'))[5]
        self.remove(H1)
        self.remove(H2)
        self.remove(H3)
        self.remove(H4)
        self.remove(H5)
        self.remove(H6)


class BT(mb.Compound):
    def __init__(self):
        super(BT,self).__init__()
        self.add(mb.load('c1ccc2nsnc2c1',smiles=True))
        carbon1 = list(self.particles_by_name('C'))[2]
        carbon2 = list(self.particles_by_name('C'))[5]
        self.add(mb.Port(anchor=carbon1,orientation=[0,1,0],separation=0.08),"p1")
        self.add(mb.Port(anchor=carbon2,orientation=[0,-1,0],separation=0.08),"p2")
        H1=list(self.particles_by_name('H'))[2]
        H2=list(self.particles_by_name('H'))[3]
        self.remove(H1)
        self.remove(H2)


class DPP(mb.Compound):
    def __init__(self):
        super(DPP, self).__init__()
        self.add(mb.load("N1C=C(C(=O)Nc2)c2C1(=O)", smiles=True))
        C1 = list(self.particles_by_name('C'))[0]
        C2 = list(self.particles_by_name('C'))[3]
        N1 = list(self.particles_by_name('N'))[0]
        N2 = list(self.particles_by_name('N'))[1]
        self.add(mb.Port(anchor=C1, orientation=[1,-3,0], separation=0.08), "p1")
        self.add(mb.Port(anchor=C2, orientation=[-1,3,0], separation=0.08), "p2")
        self.add(mb.Port(anchor=N1, orientation=[10, 1, 1], separation=0.08), "p3")
        self.add(mb.Port(anchor=N2, orientation=[-10, -1, 1], separation=0.08), "p4")
        self.remove(self.particles_by_name('H')) #doing it this way because theres no H's left on DPP


class FBT(mb.Compound):
    def __init__(self):
        super(FBT, self).__init__()
        self.add(mb.load('c1(F)ccc2nsnc2c1',smiles=True))
        carbon1 = list(self.particles_by_name('C'))[2]
        carbon2 = list(self.particles_by_name('C'))[5]
        self.add(mb.Port(anchor=carbon1, orientation=[-1, -10, 0], separation=0.08), "p1")
        self.add(mb.Port(anchor=carbon2, orientation=[1, 10, 0], separation=0.08), "p2")
        H1 = list(self.particles_by_name('H'))[1]
        H2 = list(self.particles_by_name('H'))[2]
        self.remove(H1)
        self.remove(H2)


class thiophene(mb.Compound):
    def __init__(self):
        super(thiophene, self).__init__()
        self.add(mb.load("c1ccsc1", smiles=True))
        carbon1 = list(self.particles_by_name('C'))[2]
        carbon2 = list(self.particles_by_name('C'))[3]
        self.add(mb.Port(anchor=carbon1, orientation=[-10, 1, 1], separation=0.08), "p1")
        self.add(mb.Port(anchor=carbon2, orientation=[10, -1, -1], separation=0.08), "p2")
        hydrogen1 = list(self.particles_by_name('H'))[2]
        hydrogen2 = list(self.particles_by_name('H'))[3]
        self.remove(hydrogen1)
        self.remove(hydrogen2)


class pyridine(mb.Compound):
    def __init__(self):
        super(pyridine,self).__init__()
        self.add(mb.load('c1ccncc1',smiles=True))
        carbon1 = list(self.particles_by_name('C'))[1]
        carbon2 = list(self.particles_by_name('C'))[3]
        self.add(mb.Port(anchor=carbon1, orientation=[-1,1,-1], separation=0.08), "p1")
        self.add(mb.Port(anchor=carbon2, orientation=[1,-1,1], separation=0.08), "p2")
        H1=list(self.particles_by_name('H'))[1]
        H2=list(self.particles_by_name('H'))[3]
        self.remove(H1)
        self.remove(H2)


class ene_CPDT(mb.Compound):
    def __init__(self):
        super(ene_CPDT, self).__init__()
        self.add(mb.load("C1(=C)C3=C(SC=C3)C2=C1C=CS2", smiles=True))
        carbon1 = list(self.particles_by_name('C'))[4]
        carbon2 = list(self.particles_by_name('C'))[9]
        carbon3 = list(self.particles_by_name('C'))[1]
        self.add(mb.Port(anchor=carbon1, orientation=[10, 1, 1], separation=0.08), "p1")
        self.add(mb.Port(anchor=carbon2, orientation=[-10, -1, -1], separation=0.08), "p2")
        self.add(mb.Port(anchor=carbon3, orientation=[-10, -1, -1], separation=0.08), "p3")
        self.add(mb.Port(anchor=carbon3, orientation=[10, 1, 1], separation=0.08), "p4")
        hydrogen1 = list(self.particles_by_name('H'))[2]
        hydrogen2 = list(self.particles_by_name('H'))[5]
        hydrogen3 = list(self.particles_by_name('H'))[0]
        hydrogen4 = list(self.particles_by_name('H'))[1]
        self.remove(hydrogen1)
        self.remove(hydrogen2)
        self.remove(hydrogen3)
        self.remove(hydrogen4)
