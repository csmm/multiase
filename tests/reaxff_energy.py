import pickle
from ase.structure import molecule as mol
from ase import Atoms
from multiasecalc.lammps.reaxff import ReaxFF
from multiasecalc.lammps.dynamics import LAMMPSOptimizer
from multiasecalc.utils import get_datafile

from energy import Fragmentation

from ase.data.s22 import data as s22_sim_data
from ase.data.s22 import get_number_of_dimer_atoms

from ase.data import s22

import traceback

##atomization GPAW calculator
class ReaxFFSystem():
    def __init__(self,name, atoms = None, fragment_list=None,
            minimize = False, ff_file_path=get_datafile('ffield.reax')):
        self.name = name
        if atoms:
            self.system = atoms
        else:
            self.system = mol(name) #use the default g22
        self.system.center(vacuum=10.0)
        self.calc = None
        if fragment_list:
            self.fragment_list = fragment_list
        else:
            self.fragment_list = self.system.get_chemical_symbols()
            
        self.minimize = minimize
        self.ff_file_path = ff_file_path
        
    def setup_calculator(self):          
        parameters = dict(
            neighbor   = '2.0 nsq',
            )
 
        calc = ReaxFF(ff_file_path=self.ff_file_path, parameters=parameters)
        
        return calc
    
    def get_potential_energy(self):
        self.system.set_calculator(self.setup_calculator())
        try:
            if self.minimize:
                optimizer = LAMMPSOptimizer(self.system)
                optimizer.run()
            e = self.system.get_potential_energy()
            self.system._del_calculator()
        except:
            traceback.print_exc()
            print "{0} molecule not converged".format(self.name)
            e = 1000 #not converged value
        return e


def test_s22():
    fragTest = Fragmentation(s22.s22)
    
    minimize = False
    #minimize = True
    print 'Relaxed:', minimize
    
    fragTest.molecules = []
    fragTest.fragments = []
    
    for moleculeName in s22.s22:
        data = s22.data[moleculeName]
        atoms = Atoms(data['symbols'], data['positions'])
        dimer1End = data['dimer atoms'][0]
        frag1Atoms = atoms.copy()[:dimer1End]
        frag2Atoms = atoms.copy()[dimer1End:]
        fragment1 = ReaxFFSystem(moleculeName + '_f1', frag1Atoms, minimize = minimize)
        fragment2 = ReaxFFSystem(moleculeName + '_f2', frag2Atoms, minimize = minimize)
        system = ReaxFFSystem(moleculeName, atoms, minimize = minimize,
                    fragment_list = [fragment1.name, fragment2.name])
        fragTest.molecules.append(system)
        fragTest.fragments += [fragment1, fragment2]
    
    fragTest.fill_data_reference(data_type='s22')
    fragTest.run(write=True)
    

def testSingle(molecule): 
    #TEST fragmentation
    print "Test fragmentation with s22 set"
    fragment_names = [molecule+'_f1',molecule+'_f2']
    test_f = Fragmentation([molecule])
    
    minimize = False
    #minimize = True
    print 'Relaxed:', minimize
    
    sys = Atoms(s22_sim_data[molecule]['symbols'],
                s22_sim_data[molecule]['positions'])
        
    test_f.molecules = [ReaxFFSystem(molecule,atoms=sys,
                                    fragment_list=fragment_names, minimize=minimize)]

    n = s22_sim_data[molecule]['dimer atoms']

    atom_fragments = [ sys.copy()[range(n[0])], sys.copy()[range(n[0],n[0]+n[1])] ]
    test_f.fragments = [ ReaxFFSystem(name,atoms=a, minimize = minimize)
                         for name,a in zip(fragment_names,atom_fragments)]
        
    test_f.fill_data_reference(data_type='s22')
    test_f.run(write=True)

    print test_f.data

if __name__ == "__main__":
	import sys
	
	if len(sys.argv) == 2:
		testSingle(sys.argv[1])
	else:
		test_s22()
    

