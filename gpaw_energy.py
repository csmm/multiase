import pickle
from ase.structure import molecule as mol
from ase import Atoms
from gpaw import GPAW
from gpaw.occupations import FermiDirac
from gpaw.cluster import Cluster

from energy import Fragmentation

from ase.data.s22 import data as s22_sim_data
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
    test.run()
    
    error = 5.64039227514 - test.mae
    print 'Test error: {0} eV'.format(error)

    #TEST fragmentation
    print "Test fragmentation with s22 set"
    molecule_names = ['Benzene-water_complex']
    test_f = Fragmentation(molecule_names)

    for molecule in molecule_names:
        ss = Cluster(Atoms(s22_sim_data[molecule]['symbols'],
                           s22_sim_data[molecule]['positions']))
        # split the structures
        s1 = ss.find_connected(0)
        s1.name = 's1'
        s2 = ss.find_connected(-1)
        s2.name = 's2'
        assert len(ss) == len(s1) + len(s2)    


        test_f.molecules = [ GPAWSystems(name,atoms=ss,h=0.6,vacuum=2.0,xc='PBE',
                                         fragment_list=['s1','s2']) 
                             for name in molecule_names]
        test_f.fragments = [ GPAWSystems(a.name,atoms=a,h=0.4,vacuum=2.0) for a in [s1,s2]]

#    test_f.fill_data_reference(data_type='s22')
        test_f.run(write=True)

        error =  1.1427 - test_f.mae
        print 'Test error: {0} eV'.format(error)



if __name__ == "__main__":
    main()

