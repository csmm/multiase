import pickle
from ase.structure import molecule as mol
from ase import Atoms
#from ase.calculators.lammps import LAMMPS
from lammpsCalculators.reaxff import ReaxFF

from energy import Fragmentation

from ase.data.s22 import data as s22_sim_data
from ase.data.s22 import get_number_of_dimer_atoms

from ase.data import s22

import traceback

##atomization GPAW calculator
class ReaxFFSystem():
    def __init__(self,name, atoms = None,vacuum=6.0, h=0.2, fragment_list=None, minimize = False):
        self.name = name
        if atoms:
            self.system = atoms
        else:
            self.system = mol(name) #use the default g22
        self.system.center(vacuum=vacuum)
        self.h = h
        cell = self.system.get_cell()
        self.system.set_cell((cell / (4 * h)).round() * 4 * h)
        self.system.center()
        self.calc = None
        if fragment_list:
            self.fragment_list = fragment_list
        else:
            self.fragment_list = self.system.get_chemical_symbols()
            
        self.minimize = minimize

    def setup_calculator(self):
        """
        parameters = dict(
            units      = 'real',
            atom_style = 'charge',
            pair_style = 'reax',
            pair_coeff = ['* * ffield.reax 1 2 3 4 5'],
            mass       = ['1 12.0107', '2 1.00794', '3 15.9994', '4 14.0067', '5 32.065'],
            neighbor   = '2.0 nsq'
            )
            
        if self.minimize:
            parameters['minimize'] = '1.0e-4 1.0e-6 1000 10000'
            
            
        calc = LAMMPS(
            specorder = ('C', 'H', 'O', 'N', 'S'),
            parameters = parameters,
            files = ['ffield.reax'],
            keep_tmp_files = True
            )
            """
            
        parameters = dict(
            neighbor   = '2.0 nsq',
            )
        
        if self.minimize:
            parameters['minimize'] = '1.0e-4 1.0e-6 1000 10000'
            
        calc = ReaxFF(specorder = ['C', 'H', 'O', 'N', 'S'], keep_tmp_files=True, parameters=parameters)
        return calc
    
    def get_potential_energy(self):
        self.system.set_calculator(self.setup_calculator())
        try:
            e = self.system.get_potential_energy()
            self.system._del_calculator()
        except:
            traceback.print_exc()
            print "{0} molecule not converged".format(self.name)
            e = 1000 #not converged value
        return e
    
#    def get_chemical_symbols(self):
#        return self.system.get_chemical_symbols()


def test_s22():
    fragTest = Fragmentation(s22.s22)
    
    minimize = True
    
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
    minimize = True
    
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

    print 'Relaxed:', minimize
    print test_f.data

if __name__ == "__main__":
	import sys
	
	if len(sys.argv) == 2:
		testSingle(sys.argv[1])
	else:
		test_s22()
    

