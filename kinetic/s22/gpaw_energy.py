import pickle
from ase.structure import molecule as mol
from ase import Atoms
from gpaw import GPAW
from gpaw.occupations import FermiDirac
from gpaw.cluster import Cluster
from ase.optimize import QuasiNewton

from energy import Fragmentation
from ase.data.s22 import s22 as s22_molecules
from ase.data.s22 import data as s22_sim_data
from ase.data.s22 import get_number_of_dimer_atoms
##atomization GPAW calculator
class GPAWSystems():
    def __init__(self,name, atoms = None,vacuum=6.0, h=0.17, xc='PBE',
                 setups='paw', mode='fd',fragment_list=None):
        self.name = name
        if atoms:
            self.system = atoms
        else:
            self.system = mol(name) #use the default g22
        self.system.center(vacuum=vacuum)
        cell = self.system.get_cell()
        self.h = h
        self.system.set_cell((cell / (4 * h)).round() * 4 * h)
        self.system.center()
        self.xc = xc
        self.setups = setups
        self.mode = mode
        self.calc = None
        if fragment_list:
            self.fragment_list = fragment_list
        else:
            self.fragment_list = self.system.get_chemical_symbols()

    def setup_calculator(self):
        hund = (len(self.system) == 1)
        self.system.set_calculator( GPAW(xc=self.xc,
                                         h=self.h,
                                         hund=hund,
                                         setups=self.setups,
                                         txt=self.name+'.txt',
                                         mode=self.mode,
                                         spinpol=True,
                                         nbands=-10,
                                         occupations=FermiDirac(0,
                                                                fixmagmom=True)
                                         )
                                    )
    
    def get_potential_energy(self):
        self.setup_calculator()
        try:
            relax = QuasiNewton(self.system)
            relax.run(fmax=0.05)
            e = self.system.get_potential_energy()
        except:
            print "{0} molecule not converged".format(self.name)
            e = 1000 #not converged value
        self.system._del_calculator()
        return e
    
#    def get_chemical_symbols(self):
#        return self.system.get_chemical_symbols()
def fragment_names(molecule):
    return [molecule+'_f1',molecule+'_f2']
def fragment_n(molecule,i):
    n = s22_sim_data[molecule]['dimer atoms']
    if i==1:
        r = range(n[0])
    elif i==2:
        r = range(n[0],n[0]+n[1])
    return r

def run_s22(reload=False):
    h = 0.17
    vacuum = 6.0
    #molecules = ['Ammonia_dimer']
    molecules = s22_molecules
    test_f = Fragmentation(molecules)
    
    sys = [Atoms(s22_sim_data[molecule]['symbols'],
                 s22_sim_data[molecule]['positions'])
           for molecule in molecules]
    
    test_f.molecules = [GPAWSystems(molecule,atoms=atoms,
                                    h=h, vacuum=vacuum,
                                    fragment_list=fragment_names(molecule))
                        for molecule,atoms in zip(molecules,sys)]

    atom_fragments_1 = [ system.system.copy()[fragment_n(system.name,1)] 
                         for system in test_f.molecules]
    atom_fragments_2 = [ system.system.copy()[fragment_n(system.name,2)]
                         for system in test_f.molecules]
    
    test_f.fragments = [ GPAWSystems(name+'_f1',atoms=a,h=h, vacuum=vacuum)
                         for name,a in zip(molecules,atom_fragments_1)]
    
    test_f.fragments += [ GPAWSystems(name+'_f2',atoms=a,h=h, vacuum=vacuum)
                          for name,a in zip(molecules,atom_fragments_2)]
    
    test_f.fill_data_reference(data_type='s22')
    if reload:
        test_f.load_data('data.pkl')
        test_f.calculate_error()
        test_f.calculate_mae()
        print test_f.mae
    else:
        test_f.run(write=True)



run_s22()
