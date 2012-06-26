import pickle
from ase.structure import molecule as mol
from ase import Atoms
from csmmcalc.lammps.charmm import CHARMM
from csmmcalc.lammps.dynamics import LAMMPSOptimizer
from csmmcalc.utils import get_datafile

from energy import Fragmentation

from ase.data.s22 import data as s22_sim_data
from ase.data.s22 import get_number_of_dimer_atoms

from ase.data import s22

import numpy as np
import traceback

##atomization GPAW calculator
class CHARMMSystem():
    def __init__(self,name, atoms = None, fragment_list=None, minimize = False):
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

    def setup_calculator(self):
        parameters = dict(
            neighbor   = '2.0 nsq', # bin mode seems to fail with dimers
            )
            
        calc = CHARMM(ff_file_path=get_datafile('par_all36_cgenff.prm'), parameters=parameters)
        calc._custom_thermo_args += ['evdwl', 'ecoul', 'emol', 'ebond']
        return calc
    
    def get_potential_energy(self):
        self.system.set_calculator(self.setup_calculator())
        try:
            if self.minimize:
                optimizer = LAMMPSOptimizer(self.system)
                optimizer.run(min_style='cg')
            e = self.system.get_potential_energy()
            self.system._del_calculator()
        except:
            traceback.print_exc()
            print "{0} molecule not converged".format(self.name)
            e = np.nan #not converged value
        return e

def calculate_charges(atoms):
    from csmmcalc.lammps.reaxff import ReaxFF
    atoms.center(vacuum=1)
    atoms.calc = ReaxFF(ff_file_path=get_datafile('ffield.reax'), parameters = dict(neighbor='2.0 nsq'))
    atoms.get_potential_energy()

def test_s22():
    fragTest = Fragmentation(s22.s22)
    
    minimize = False
    #minimize = True
    print 'Relaxed:', minimize
    
    fragTest.molecules = []
    fragTest.fragments = []
    
    testSet = [
        'Ammonia_dimer',
        'Water_dimer',
        'Formic_acid_dimer',
        'Formamide_dimer',
        'Uracil_dimer_h-bonded',
        #'2-pyridoxine_2-aminopyridine_complex',
        'Adenine-thymine_Watson-Crick_complex',
        'Methane_dimer',
        'Ethene_dimer',
        'Benzene-methane_complex',
        'Benzene_dimer_parallel_displaced',
        'Pyrazine_dimer',
        'Uracil_dimer_stack',
        'Indole-benzene_complex_stack',
        'Adenine-thymine_complex_stack',
        #'Ethene-ethyne_complex',
        'Benzene-water_complex',
        'Benzene-ammonia_complex',
        #'Benzene-HCN_complex',
        'Benzene_dimer_T-shaped',
        'Indole-benzene_T-shape_complex',
        'Phenol_dimer'
        ]
    
    for moleculeName in testSet:
        atoms = s22.create_s22_system(moleculeName)
        dimer1End = s22.data[moleculeName]['dimer atoms'][0]
        frag1Atoms = atoms.copy()[:dimer1End]
        frag2Atoms = atoms.copy()[dimer1End:]
        calculate_charges(frag1Atoms)
        calculate_charges(frag2Atoms)
        atoms.set_charges(np.append(frag1Atoms.get_charges(), frag2Atoms.get_charges()))
        fragment1 = CHARMMSystem(moleculeName + '_f1', frag1Atoms, minimize = minimize)
        fragment2 = CHARMMSystem(moleculeName + '_f2', frag2Atoms, minimize = minimize)
        system = CHARMMSystem(moleculeName, atoms, minimize = minimize,
                    fragment_list = [fragment1.name, fragment2.name])
        fragTest.molecules.append(system)
        fragTest.fragments += [fragment1, fragment2]
    
    fragTest.fill_data_reference(data_type='s22')
    fragTest.run(write=True)
    

def testSingle(moleculeName = "Methane_dimer"): 
    #TEST fragmentation
    print "Test fragmentation with s22 set"
    
    minimize = False
    minimize = True
    print 'Relaxed:', minimize
    
    fragTest = Fragmentation([moleculeName])
    
    #atoms = s22.create_s22_system(moleculeName)
    data = s22.data[moleculeName]
    atoms = Atoms(data['symbols'], data['positions'])
    dimer1End = s22.data[moleculeName]['dimer atoms'][0]
    frag1Atoms = atoms.copy()[:dimer1End]
    frag2Atoms = atoms.copy()[dimer1End:]
    
    calculate_charges(frag1Atoms)
    calculate_charges(frag2Atoms)
    atoms.set_charges(np.append(frag1Atoms.get_charges(), frag2Atoms.get_charges()))
    
    fragment1 = CHARMMSystem(moleculeName + '_f1', frag1Atoms, minimize = minimize)
    fragment2 = CHARMMSystem(moleculeName + '_f2', frag2Atoms, minimize = minimize)
    system = CHARMMSystem(moleculeName, atoms, minimize = minimize,
                fragment_list = [fragment1.name, fragment2.name])
    fragTest.molecules = [system]
    fragTest.fragments = [fragment1, fragment2]
        
    fragTest.fill_data_reference(data_type='s22')
    fragTest.run(write=True)
    
    import ase.io    
    ase.io.write('molecule.xyz', atoms)

    print fragTest.data


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2:
        testSingle(sys.argv[1])
    else:
        test_s22()

