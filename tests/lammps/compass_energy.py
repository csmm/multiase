import pickle
from ase.structure import molecule as mol
from ase import Atoms
from multiasecalc.lammps.compass import COMPASS
from multiasecalc.lammps.dynamics import LAMMPSOptimizer
from multiasecalc.utils import get_datafile

from energy import Fragmentation

from ase.data.s22 import data as s22_sim_data
from ase.data.s22 import get_number_of_dimer_atoms

from ase.data import s22

import numpy as np
import traceback

##atomization GPAW calculator
class COMPASSSystem():
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
        parameters = dict(
            neighbor   = '2.0 nsq', # bin mode seems to fail with dimers
            )
            
        calc = COMPASS(ff_file_path=get_datafile('compass.frc'), parameters=parameters, debug=True)
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
            e = np.nan #not converged value
        return e


def test_s22():
    fragTest = Fragmentation(s22.s22)
    
    minimize = False
    minimize = True
    print 'Relaxed:', minimize
    
    fragTest.molecules = []
    fragTest.fragments = []
    
    testSet = [
		#'Ammonia_dimer',
		#'Water_dimer',
		#'Formic_acid_dimer',
		#'Formamide_dimer',
		#'Uracil_dimer_h-bonded',
		#'2-pyridoxine_2-aminopyridine_complex',
		#'Adenine-thymine_Watson-Crick_complex',
		'Methane_dimer',
		#'Ethene_dimer',
		'Benzene-methane_complex',
		#'Benzene_dimer_parallel_displaced',
		#'Pyrazine_dimer',
		#'Uracil_dimer_stack',
		#'Indole-benzene_complex_stack',
		#'Adenine-thymine_complex_stack',
		#'Ethene-ethyne_complex',
		#'Benzene-water_complex',
		#'Benzene-ammonia_complex',
		#'Benzene-HCN_complex',
		#'Benzene_dimer_T-shaped',
		#'Indole-benzene_T-shape_complex',
		'Phenol_dimer'
		]
    
    for moleculeName in testSet:
        data = s22.data[moleculeName]
        atoms = Atoms(data['symbols'], data['positions'])
        dimer1End = data['dimer atoms'][0]
        frag1Atoms = atoms.copy()[:dimer1End]
        frag2Atoms = atoms.copy()[dimer1End:]
        fragment1 = COMPASSSystem(moleculeName + '_f1', frag1Atoms, minimize = minimize)
        fragment2 = COMPASSSystem(moleculeName + '_f2', frag2Atoms, minimize = minimize)
        system = COMPASSSystem(moleculeName, atoms, minimize = minimize,
                    fragment_list = [fragment1.name, fragment2.name])
        fragTest.molecules.append(system)
        fragTest.fragments += [fragment1, fragment2]
    
    fragTest.fill_data_reference(data_type='s22')
    fragTest.run(write=True)
    

def testSingle(molecule = "Methane_dimer"): 
    #TEST fragmentation
    print "Test fragmentation with s22 set"
    fragment_names = [molecule+'_f1',molecule+'_f2']
    test_f = Fragmentation([molecule])
    
    minimize = False
    minimize = True
    print 'Relaxed:', minimize
    
    sys = Atoms(s22_sim_data[molecule]['symbols'],
                s22_sim_data[molecule]['positions'])
        
    test_f.molecules = [COMPASSSystem(molecule,atoms=sys,
                                    fragment_list=fragment_names, minimize=minimize)]

    n = s22_sim_data[molecule]['dimer atoms']

    atom_fragments = [ sys.copy()[range(n[0])], sys.copy()[range(n[0],n[0]+n[1])] ]
    test_f.fragments = [ COMPASSSystem(name,atoms=a, minimize = minimize)
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

