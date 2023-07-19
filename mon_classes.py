import mbuild as mb
import warnings 
import numpy
from frag_classes import IDT, CPDT, DPP, BT, PT, FBT, thiophene, pyridine, TPD,BDT
from r_classes import c11_bo,HD,ODD,C1BO,C3BO,C4BO,C5BO, C16, ene_HD, ene_ODD
import ipywidgets as widgets
warnings.filterwarnings('ignore')


class PIDTCPDT_C11BO(mb.Compound):
    def __init__(self, energy_minimize=False):
        super(PIDTCPDT_C11BO, self).__init__()
        idt=IDT()
        cpdt=CPDT()
        r1 = c11_bo()
        r2 = c11_bo()
        r3 = c11_bo()
        r4 = c11_bo()
        r5 = c11_bo()
        r6 = c11_bo()
        self.add([idt,cpdt,r1,r2,r3,r4,r5,r6])
        mb.force_overlap(move_this=idt,
                        to_positions=cpdt['p1'],
                        from_positions=idt['p1'])
        mb.force_overlap(move_this=r1,
                        to_positions=cpdt['p3'],
                        from_positions=r1['p1'])
        mb.force_overlap(move_this=r2,
                        to_positions=cpdt['p4'],
                        from_positions=r2['p1'])
        mb.force_overlap(move_this=r3,
                        to_positions=idt['p3'],
                        from_positions=r3['p1'])
        mb.force_overlap(move_this=r4,
                        to_positions=idt['p4'],
                        from_positions=r4['p1'])
        mb.force_overlap(move_this=r5,
                        to_positions=idt['p5'],
                        from_positions=r5['p1'])
        mb.force_overlap(move_this=r6,
                        to_positions=idt['p6'],
                        from_positions=r6['p1'])
        r1.translate([0,0,2])
        r2.translate([0,0,2])
        r3.translate([0,0,2])
        r4.translate([0,0,2])
        r5.translate([0,0,2])
        r6.translate([0,0,2])
        cpdt.translate([0,1,0])
        self.bond_indices = [13,26]
        self.orientations = [[-1,0,0],[1,0,0]]
        self.separation = 0.14
        self.replace = False
        if energy_minimize == True:
            self.energy_minimize()



class P3HT(mb.Compound):
    def __init__(self):
        super(P3HT,self).__init__()
        self.add(mb.load("CCCCCCC1=C(SC(=C1))",smiles=True))
        self.bond_indices = [24,25]
        self.orientations = [[0,0,1],[0,0,-1]]
        self.separation = 0.14
        self.replace = True 
       # if energy_minimize == True:
       #     self.energy_minimize()


class PCPDTPT_HD(mb.Compound):
    def __init__(self):
        super(PCPDTPT_HD,self).__init__()
        cpdt=CPDT()
        pt=PT()
        r1=HD()
        r2 = HD()
        self.add([cpdt,pt,r1,r2])
        mb.force_overlap(move_this=pt,
                        to_positions=cpdt['p1'],
                        from_positions=pt['p2'])
        mb.force_overlap(move_this=r1,
                        to_positions=cpdt['p3'],
                        from_positions=r1['p1'])
        mb.force_overlap(move_this=r2,
                        to_positions=cpdt['p4'],
                        from_positions=r2['p1'])
        self.bond_indices = [4,21]
        self.orientations = [[-1,0,0],[1,0,0]]
        self.separation = 0.14
        self.replace = False        


class PCPDTPT_ODD(mb.Compound):
    def __init__(self):
        super(PCPDTPT_ODD,self).__init__()
        cpdt=CPDT()
        pt=PT()
        r1=ODD()
        r2 = ODD()
        self.add([cpdt,pt,r1,r2])
        mb.force_overlap(move_this=pt,
                        to_positions=cpdt['p1'],
                        from_positions=pt['p2'])
        mb.force_overlap(move_this=r1,
                        to_positions=cpdt['p3'],
                        from_positions=r1['p1'])
        mb.force_overlap(move_this=r2,
                        to_positions=cpdt['p4'],
                        from_positions=r2['p1'])
        self.bond_indices = [4,21]
        self.orientations = [[-1,0,0],[1,0,0]]
        self.separation = 0.14
        self.replace = False        

class PDPPPyT_ODD(mb.Compound):
    def __init__(self):
        super(PDPPPyT_ODD,self).__init__()
        dpp=DPP()
        pyr1=pyridine()
        pyr2=pyridine()
        T1=thiophene()
        T2=thiophene()
        r1=ODD()
        r2 = ODD()
        self.add([dpp,T1,T2,pyr1,pyr2,r1,r2])
        mb.force_overlap(move_this=pyr1,
                        to_positions=T1['p1'],
                        from_positions=pyr1['p1'])
        mb.force_overlap(move_this=dpp,
                        to_positions=pyr1['p2'],
                        from_positions=dpp['p1'])
        mb.force_overlap(move_this=pyr2,
                        to_positions=dpp['p2'],
                        from_positions=pyr2['p2'])
        mb.force_overlap(move_this=T2,
                        to_positions=pyr2['p1'],
                        from_positions=T2['p1'])
        mb.force_overlap(move_this=r1,
                        to_positions=dpp['p3'],
                        from_positions=r1['p1'])
        mb.force_overlap(move_this=r2,
                        to_positions=dpp['p4'],
                        from_positions=r2['p1'])
        self.bond_indices = [14,21]
        self.orientations = [[1,1,1],[-1,-1,-1]]
        self.separation = 0.14
        self.replace = False 


class PCPDTFBT_C1_BO(mb.Compound):
    def __init__(self):
        super(PCPDTFBT_C1_BO,self).__init__()
        cpdt1=CPDT()
        cpdt2=CPDT()
        fbt1=FBT()
        fbt2=FBT()
        r1=C1BO()
        r2=C1BO()
        r3=C1BO()
        r4=C1BO()
        self.add([cpdt1,cpdt2,fbt1,fbt2,r1,r2,r3,r4])
        mb.force_overlap(move_this=fbt1,
                        to_positions=cpdt1['p2'],
                        from_positions=fbt1['p1'])
        mb.force_overlap(move_this=cpdt2,
                        to_positions=fbt1['p2'],
                        from_positions=cpdt2['p1'])
        mb.force_overlap(move_this=fbt2,
                        to_positions=cpdt2['p2'],
                        from_positions=fbt2['p2'])
        mb.force_overlap(move_this=r1,
                        to_positions=cpdt1['p3'],
                        from_positions=r1['p1'])
        mb.force_overlap(move_this=r2,
                        to_positions=cpdt1['p4'],
                        from_positions=r2['p1'])
        mb.force_overlap(move_this=r3,
                        to_positions=cpdt2['p3'],
                        from_positions=r3['p1'])
        mb.force_overlap(move_this=r4,
                        to_positions=cpdt2['p4'],
                        from_positions=r4['p1'])
        self.bond_indices = [9,40]
        self.orientations = [[-1,0,0],[1,0,0]]
        self.separation = 0.14
        self.replace = False


class PCPDTFBT_C3_BO(mb.Compound):
    def __init__(self):
        super(PCPDTFBT_C3_BO,self).__init__()
        cpdt1=CPDT()
        cpdt2=CPDT()
        fbt1=FBT()
        fbt2=FBT()
        r1=C3BO()
        r2=C3BO()
        r3=C3BO()
        r4=C3BO()
        self.add([cpdt1,cpdt2,fbt1,fbt2,r1,r2,r3,r4])
        mb.force_overlap(move_this=fbt1,
                        to_positions=cpdt1['p2'],
                        from_positions=fbt1['p1'])
        mb.force_overlap(move_this=cpdt2,
                        to_positions=fbt1['p2'],
                        from_positions=cpdt2['p1'])
        mb.force_overlap(move_this=fbt2,
                        to_positions=cpdt2['p2'],
                        from_positions=fbt2['p2'])
        mb.force_overlap(move_this=r1,
                        to_positions=cpdt1['p3'],
                        from_positions=r1['p1'])
        mb.force_overlap(move_this=r2,
                        to_positions=cpdt1['p4'],
                        from_positions=r2['p1'])
        mb.force_overlap(move_this=r3,
                        to_positions=cpdt2['p3'],
                        from_positions=r3['p1'])
        mb.force_overlap(move_this=r4,
                        to_positions=cpdt2['p4'],
                        from_positions=r4['p1'])
        r1.rotate(theta=1.0,around=[0,1,0])
        r2.rotate(theta=1.0,around=[0,1,0])
        r3.rotate(theta=1.0,around=[0,1,0])
        r4.rotate(theta=1.0,around=[0,1,0])
        self.bond_indices = [9,40]
        self.orientations = [[-1,0,0],[1,0,0]]
        self.separation = 0.14
        self.replace = False

class PCPDTFBT_C4_BO(mb.Compound):
    def __init__(self):
        super(PCPDTFBT_C4_BO,self).__init__()
        cpdt1=CPDT()
        cpdt2=CPDT()
        fbt1=FBT()
        fbt2=FBT()
        r1=C4BO()
        r2=C4BO()
        r3=C4BO()
        r4=C4BO()
        self.add([cpdt1,cpdt2,fbt1,fbt2,r1,r2,r3,r4])
        mb.force_overlap(move_this=fbt1,
                        to_positions=cpdt1['p2'],
                        from_positions=fbt1['p1'])
        mb.force_overlap(move_this=cpdt2,
                        to_positions=fbt1['p2'],
                        from_positions=cpdt2['p1'])
        mb.force_overlap(move_this=fbt2,
                        to_positions=cpdt2['p2'],
                        from_positions=fbt2['p2'])
        mb.force_overlap(move_this=r1,
                        to_positions=cpdt1['p3'],
                        from_positions=r1['p1'])
        mb.force_overlap(move_this=r2,
                        to_positions=cpdt1['p4'],
                        from_positions=r2['p1'])
        mb.force_overlap(move_this=r3,
                        to_positions=cpdt2['p3'],
                        from_positions=r3['p1'])
        mb.force_overlap(move_this=r4,
                        to_positions=cpdt2['p4'],
                        from_positions=r4['p1'])
        r1.rotate(theta=1.0,around=[0,1,0])
        r2.rotate(theta=1.0,around=[0,1,0])
        r3.rotate(theta=1.0,around=[0,1,0])
        r4.rotate(theta=1.0,around=[0,1,0])
        self.bond_indices = [9,40]
        self.orientations = [[-1,0,0],[1,0,0]]
        self.separation = 0.14
        self.replace = False


class PCPDTFBT_C5_BO(mb.Compound):
    def __init__(self):
        super(PCPDTFBT_C5_BO,self).__init__()
        cpdt1=CPDT()
        cpdt2=CPDT()
        fbt1=FBT()
        fbt2=FBT()
        r1=C5BO()
        r2=C5BO()
        r3=C5BO()
        r4=C5BO()
        self.add([cpdt1,cpdt2,fbt1,fbt2,r1,r2,r3,r4])
        mb.force_overlap(move_this=fbt1,
                        to_positions=cpdt1['p2'],
                        from_positions=fbt1['p1'])
        mb.force_overlap(move_this=cpdt2,
                        to_positions=fbt1['p2'],
                        from_positions=cpdt2['p1'])
        mb.force_overlap(move_this=fbt2,
                        to_positions=cpdt2['p2'],
                        from_positions=fbt2['p2'])
        mb.force_overlap(move_this=r1,
                        to_positions=cpdt1['p3'],
                        from_positions=r1['p1'])
        mb.force_overlap(move_this=r2,
                        to_positions=cpdt1['p4'],
                        from_positions=r2['p1'])
        mb.force_overlap(move_this=r3,
                        to_positions=cpdt2['p3'],
                        from_positions=r3['p1'])
        mb.force_overlap(move_this=r4,
                        to_positions=cpdt2['p4'],
                        from_positions=r4['p1'])
        r1.rotate(theta=1.0,around=[0,1,0])
        r2.rotate(theta=1.0,around=[0,1,0])
        r3.rotate(theta=1.0,around=[0,1,0])
        r4.rotate(theta=1.0,around=[0,1,0])
        self.bond_indices = [9,40]
        self.orientations = [[-1,0,0],[1,0,0]]
        self.separation = 0.14
        self.replace = False


class PCPDTFBT_C11_BO(mb.Compound):
    def __init__(self):
        super(PCPDTFBT_C11_BO,self).__init__()
        cpdt1=CPDT()
        cpdt2=CPDT()
        fbt1=FBT()
        fbt2=FBT()
        r1=c11_bo()
        r2=c11_bo()
        r3=c11_bo()
        r4=c11_bo()
        self.add([cpdt1,cpdt2,fbt1,fbt2,r1,r2,r3,r4])
        mb.force_overlap(move_this=fbt1,
                        to_positions=cpdt1['p2'],
                        from_positions=fbt1['p1'])
        mb.force_overlap(move_this=cpdt2,
                        to_positions=fbt1['p2'],
                        from_positions=cpdt2['p1'])
        mb.force_overlap(move_this=fbt2,
                        to_positions=cpdt2['p2'],
                        from_positions=fbt2['p2'])
        mb.force_overlap(move_this=r1,
                        to_positions=cpdt1['p3'],
                        from_positions=r1['p1'])
        mb.force_overlap(move_this=r2,
                        to_positions=cpdt1['p4'],
                        from_positions=r2['p1'])
        mb.force_overlap(move_this=r3,
                        to_positions=cpdt2['p3'],
                        from_positions=r3['p1'])
        mb.force_overlap(move_this=r4,
                        to_positions=cpdt2['p4'],
                        from_positions=r4['p1'])
        r1.rotate(theta=1.0,around=[0,1,0])
        r2.rotate(theta=1.0,around=[0,1,0])
        r3.rotate(theta=1.0,around=[0,1,0])
        r4.rotate(theta=1.0,around=[0,1,0])
        self.bond_indices = [9,40]
        self.orientations = [[-1,0,0],[1,0,0]]
        self.separation = 0.14
        self.replace = False



class PCPDTPT_nC16(mb.Compound):
    def __init__(self):
        super(PCPDTPT_nC16,self).__init__()
        cpdt = CPDT()
        pt = PT()
        r1 = C16()
        r2 = C16()
        self.add([cpdt,pt,r1,r2])
        mb.force_overlap(move_this=pt,
                         to_positions=cpdt['p1'],
                         from_positions=pt['p1'])
        mb.force_overlap(move_this=r1,
                        to_positions=cpdt['p3'],
                        from_positions=r1['p1'])
        mb.force_overlap(move_this=r2,
                        to_positions=cpdt['p4'],
                        from_positions=r2['p1'])
        self.bond_indices = [4,15]
        self.orientations = [[1,0,0],[-1,0,0]]
        self.separation = 0.14
        self.replace = False


class PIDTBT_nC16(mb.Compound):
    def __init__(self):
        super(PIDTBT_nC16,self).__init__()
        idt = IDT()
        bt = BT()
        r1 = C16()
        r2 = C16()
        r3 = C16()
        r4 = C16()
        self.add([idt,bt,r1,r2,r3,r4])
        mb.force_overlap(move_this=bt,
                         to_positions=idt['p1'],
                         from_positions=bt['p1'])
        mb.force_overlap(move_this=r1,
                        to_positions=idt['p3'],
                        from_positions=r1['p1'])
        mb.force_overlap(move_this=r2,
                        to_positions=idt['p4'],
                        from_positions=r2['p1'])
        mb.force_overlap(move_this=r3,
                        to_positions=idt['p5'],
                        from_positions=r3['p1'])
        mb.force_overlap(move_this=r4,
                        to_positions=idt['p6'],
                        from_positions=r4['p1'])
        self.bond_indices = [13,30]
        self.orientations = [[1,0,0],[-1,0,0]]
        self.separation = 0.14
        self.replace = False


class PIDTFBT_C11_BO(mb.Compound):
    def __init__(self):
        super(PIDTFBT_C11_BO,self).__init__()
        idt1 = IDT()
        fbt1 = FBT()
        idt2 = IDT()
        fbt2 = FBT()
        r1 = c11_bo()
        r2 = c11_bo()
        r3 = c11_bo()
        r4 = c11_bo()
        r5 = c11_bo()
        r6 = c11_bo()
        r7 = c11_bo()
        r8 = c11_bo()
        self.add([idt1,fbt1,idt2,fbt2,r1,r2,r3,r4,r5,r6,r7,r8])
        mb.force_overlap(move_this=fbt1,
                         to_positions=idt1['p1'],
                         from_positions=fbt1['p1'])
        mb.force_overlap(move_this=idt2,
                         to_positions=fbt1['p2'],
                         from_positions=idt2['p1'])
        mb.force_overlap(move_this=fbt2,
                         to_positions=idt2['p2'],
                         from_positions=fbt2['p2'])
        mb.force_overlap(move_this=r1,
                         to_positions=idt1['p3'],
                         from_positions=r1['p1'])
        mb.force_overlap(move_this=r2,
                         to_positions=idt1['p4'],
                         from_positions=r2['p1'])
        mb.force_overlap(move_this=r3,
                         to_positions=idt1['p5'],
                         from_positions=r3['p1'])
        mb.force_overlap(move_this=r4,
                         to_positions=idt1['p6'],
                         from_positions=r4['p1'])
        mb.force_overlap(move_this=r5,
                         to_positions=idt2['p3'],
                         from_positions=r5['p1'])
        mb.force_overlap(move_this=r6,
                         to_positions=idt2['p4'],
                         from_positions=r6['p1'])
        mb.force_overlap(move_this=r7,
                         to_positions=idt2['p5'],
                         from_positions=r7['p1'])
        mb.force_overlap(move_this=r8,
                         to_positions=idt2['p6'],
                         from_positions=r8['p1'])
        self.bond_indices = [13,58]
        self.orientations = [[1,0,0],[-1,0,0]]
        self.separation = 0.14
        self.replace = False



class BDT_TPD(mb.Compound):
    def __init__(self):
        super(BDT_TPD,self).__init__()
        tpd = TPD()
        bdt = BDT()
        self.add([tpd,bdt])
        mb.force_overlap(move_this=bdt,
                        to_positions=tpd['p2'],
                        from_positions=bdt['p1'])
        self.bond_indices = [3,49]
        self.orientations = [[-1,0,0],[1,0,0]]
        self.separation = 0.14
        self.replace = False


class perylene(mb.Compound):
    def __init__(self):
        super(perylene,self).__init__()
        self.add(mb.load('c1cc2cccc3c4cccc5cccc(c(c1)c23)c45',smiles=True))



class PCPDTPT_eneODD(mb.Compound):
    def __init__(self):
        super(PCPDTPT_eneODD,self).__init__()
        cpdt = CPDT()
        pt = PT()
        r = ene_ODD()
        self.add([cpdt,pt,r])
        mb.force_overlap(move_this=pt,
                        to_positions=cpdt['p2'],
                        from_positions=pt['p2'])
        mb.force_overlap(move_this=r,
                        to_positions=cpdt['p3'],
                        from_positions=r['p1'])
        mb.force_overlap(move_this=r,
                        to_positions=cpdt['p4'],
                        from_positions=r['p2'])
        self.bond_indices = [9,21]
        self.orientations = [[-1,0,0],[1,0,0]]
        self.separation = 0.14
        self.replace = False
        r.translate([1,1,1])

        
class PCPDT_PT_eneHD(mb.Compound):
    def __init__(self):
        super(PCPDT_PT_eneHD,self).__init__()
        cpdt = CPDT()
        pt = PT()
        r = ene_HD()
        self.add([cpdt,pt,r])
        mb.force_overlap(move_this=pt,
                        to_positions=cpdt['p2'],
                        from_positions=pt['p2'])
        mb.force_overlap(move_this=r,
                        to_positions=cpdt['p3'],
                        from_positions=r['p1'])
        mb.force_overlap(move_this=r,
                        to_positions=cpdt['p4'],
                        from_positions=r['p2'])
        self.bond_indices = [9,21]
        self.orientations = [[-1,0,0],[1,0,0]]
        self.separation = 0.14
        self.replace = False
