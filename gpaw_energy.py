import pickle
from ase.structure import molecule as mol
from ase import Atoms
from gpaw import GPAW
from gpaw.occupations import FermiDirac
from gpaw.cluster import Cluster

from energy import Fragmentation

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
    
#    def get_chemical_symbols(self):
#        return self.system.get_chemical_symbols()



def main():
    #TEST atomization
    print "Test atomization"
    molecule_names = ['CH4']
    test = Fragmentation(molecule_names)
    test.molecules = [GPAWSystems(name,h=0.3,vacuum=2.0,xc='PBE')
                      for name in molecule_names]
    
    atom_names = []
    for gs in test.molecules:
        atom_names.extend(gs.system.get_chemical_symbols())
  
    test.fragments = [GPAWSystems(name,h=0.3,vacuum=2.0,xc='PBE') 
                      for name in set(atom_names)]

    test.fill_data_reference(data_type='g22')
    #test.run()
    
    error = 5.64039227514 - test.mae
    print 'Test error: {0} eV'.format(error)

    #TEST fragmentation
    print "Test fragmentation with s22 set"
    molecule = 'Ammonia_dimer'
    fragment_names = [molecule+'_f1',molecule+'f2']
    test_f = Fragmentation([molecule])
    
    sys = Atoms(s22_sim_data[molecule]['symbols'],
                s22_sim_data[molecule]['positions'])
    
    test_f.molecules = [GPAWSystems(molecule,atoms=sys,
                                    h=0.2, vacuum=3.0,
                                    fragment_list=fragment_names)]

    n = s22_sim_data[molecule]['dimer atoms']

    atom_fragments = [ sys.copy()[range(n[0])], sys.copy()[range(n[0],n[0]+n[1])] ]
    test_f.fragments = [ GPAWSystems(name,atoms=a,h=0.2, vacuum=3.0)
                         for name,a in zip(fragment_names,atom_fragments)]
        
    test_f.fill_data_reference(data_type='s22')
    test_f.run(write=True)
        
    error =  1.1427 - test_f.mae
    print 'Test error: {0} eV'.format(error)



if __name__ == "__main__":
    main()

