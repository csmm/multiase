import pickle
from ase.structure import molecule as mol
from gpaw import GPAW
from gpaw.occupations import FermiDirac
#from ase.data.molecules import molecule, atoms as g2_atoms, g1
from ase.data.g2_1 import atom_names as g22_atom_names 
from ase.data.g2_1 import molecule_names as g22_molecule_names 
from ase.data.g2_1 import data as g22_data
import g2_1 as short_molecules_exp_data
from energy import Atomization


##atomization GPAW calculator
class GPAWSystems():
    def __init__(self,name, vacuum=6.0, h=0.17, xc='LDA',
                 setups='paw', mode='fd'):
        self.name = name
        self.system = mol(name)
        self.system.center(vacuum=vacuum)
        cell = self.system.get_cell()
        self.h = h
        self.system.set_cell((cell / (4 * h)).round() * 4 * h)
        self.system.center()
        self.xc = xc
        self.setups = setups
        self.mode = mode
        self.calc = None
    def setup_calculator(self):
        hund = (len(self.system) == 1)
        calc = GPAW(xc=self.xc,
                    h=self.h,
                    hund=hund,
                    setups=self.setups,
                    txt=self.name+'.txt',
                    mode=self.mode,
                    spinpol=True,
                    occupations=FermiDirac(0,fixmagmom=True)
                    )
        return calc
    
    def get_potential_energy(self):
        self.system.set_calculator(self.setup_calculator())
        try:
            e = self.system.get_potential_energy()
            self.system._del_calculator()
        except:
            print "{0} molecule not converged".format(self.name)
            e = 1000 #not converged value
        return e
    
    def get_chemical_symbols(self):
        return self.system.get_chemical_symbols()



def main():
    test = Atomization()
    #calling g22 or s22
    #molecule_names = g22_molecule_names # ['H2,...]
    #atom_names = g22_atom_names #['H',...]
    #test.exp_data = g22_data
    #molecule_names = ['CH3','CH4','OH','H2O','C2H2','C2H4','C2H6','CO']

    molecule_names = ['CH4']
    atom_names = ['H','C']
    test.exp_data = short_molecules_exp_data

    test.add_systems(molecule_names)
    test.atoms = [GPAWSystems(name,h=0.3,vacuum=2.0,xc='PBE') 
                  for name in atom_names]
    test.molecules = [GPAWSystems(name,h=0.3,vacuum=2.0,xc='PBE') 
                      for name in molecule_names]
    test.fill_data_experimental(factor_to_eV=0.0433641146392)

    test.get_all_energies()
    test.write_energies('energies.pkl')
    test.fill_data_simulated()
    test.calculate_error()
    test.calculate_mae()
    test.write_data('data.pkl')
    
    error = 5.64039227514 - test.mae
    print 'Test error: {0} eV'.format(error)

if __name__ == "__main__":
    main()

