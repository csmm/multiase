import pickle
from ase.structure import molecule as mol
from gpaw import GPAW
from gpaw.occupations import FermiDirac
import g2_1 as g22_exp_data
import s22 as s22_sim_data
#from ase.data import g2_1 as exp_data

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
    print "Test atomization"
    test = Atomization()
    molecule_names = ['CH4']
    atom_names = ['H','C']
    test.ref_data = g22_exp_data

    test.add_systems(molecule_names)
    test.atoms = [GPAWSystems(name,h=0.3,vacuum=2.0,xc='PBE') 
                  for name in atom_names]
    test.molecules = [GPAWSystems(name,h=0.3,vacuum=2.0,xc='PBE') 
                      for name in molecule_names]
    test.fill_data_reference(data_type='g22')

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

