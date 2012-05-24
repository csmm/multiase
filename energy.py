import pickle
import g2_1 as g22_exp_data
from ase.data.s22 import data as s22_sim_data
from ase.data.s22 import get_interaction_energy_s22

class ScalarTest(object):
    def __init__(self):
        self.data = {}
        self.mae = 0.0
    def add_systems(self,names): 
        for name in names:
            self.data[name] = {'ref':0.0,'sim':0.0,'err':0.0}
    def add_data(self, name, data, data_type):
        """ Data types are 'ref' or 'sim' """
        self.data[name][data_type] = data
    def calculate_error(self):
        for name, values in self.data.iteritems():
            values['err'] = values['ref']-values['sim']
    def calculate_mae(self):
        for name, values in self.data.iteritems():
            self.mae += abs(values['ref']-values['sim'])
        self.mae /= len(self.data)
        #print 'Mean absolute error is {0} eV'.format(self.mae)
    def write_data(self,filename):
        f = open(filename,'wb')
        pickle.dump(self.data, f)
        f.close()
        f = open('data.txt','w')
        for k,v in self.data.iteritems():
            out = str(k)+'\t'+str(v)+'\n'
            f.write(out)
        f.write('mae: {0}'.format(self.mae))
    def load_data(self,filename):
        f = open(filename,'rb')
        self.data = pickle.load(f)
        f.close()

class Fragmentation(ScalarTest):
    def __init__(self,molecule_names):
        self.molecules = None #List of molecules
        self.fragments = None #Atom object ase-like
        self.calc = None  #Calculator gpaw-like
        self.energies = {}
        self.data = {}
        self.mae = 0.0
        self.data_ref = None
        self.add_systems(molecule_names)

    def get_all_energies(self):
        for fragment in set(self.fragments):
            e = fragment.get_potential_energy()
            self.energies[fragment.name] = e
        for system in self.molecules:
            e = system.get_potential_energy()
            self.energies[system.name] = e
    
    def fill_data_simulated(self):
        for system in self.molecules:
            ae = -1.0 * self.energies[system.name]
            for fragment in system.fragment_list:
                ae += self.energies[fragment]
            self.add_data(system.name,ae,'sim')

    def fill_data_reference(self, data_type='g22'):
        """ Energies in eV"""
        for system in self.molecules:
            if data_type == 'g22':
                to_eV=0.0433641146392
                e = g22_exp_data.get_atomization_energy(system.name) * to_eV
            elif data_type == 's22':
                e = -get_interaction_energy_s22(system.name)
            else:
                print "Data type not implemented"
            self.add_data(system.name,e,'ref')
 
    def write_energies(self,filename):
        f = open(filename,'wb')
        pickle.dump(self.energies, f)
        f.close()
        f = open('energies.txt','w')
        for k,v in self.energies.iteritems():
            out = str(k)+'\t'+str(v)+'\n'
            f.write(out)

    def load_energies(self,filename):
        f = open(filename,'rb')
        self.energies = pickle.load(f)
        f.close()

    def run(self,write=True):
        self.get_all_energies()
        self.fill_data_simulated()
        self.calculate_error()
        self.calculate_mae()
        if write:
            self.write_energies('energies.pkl')
            self.write_data('data.pkl')

        
class SimpleAtoms():        
    def __init__(self,name,fragment_list=None):
        self.name = name
        self.fragment_list = fragment_list
        return
    def get_potential_energy(self):
        return 1.0


def main():
    print "Test atomization using exp data in g22 set"
    molecule_names = ['CH4']
    fragment_names = ['H','H','H','H','C']
    test = Fragmentation(molecule_names)

    test.fragments = [ SimpleAtoms(name) for name in set(fragment_names) ]
    test.molecules = [ SimpleAtoms(name,fragment_list=fragment_names) for name in molecule_names ]

    test.fill_data_reference(data_type='g22')
    test.run()

    error = 14.220659970104172 - test.mae
    print 'Test error: {0} eV'.format(error)

    #to load
    print "Test load data"
    test_d = ScalarTest()
    test_d.load_data('data.pkl')
    test_d.calculate_mae()
    error = 14.220659970104172 - test_d.mae
    print 'Test error: {0} eV'.format(error)

    print "Test fragmentation with s22 set"
    molecule_names = ['Benzene-water_complex']
    fragment_names = ['Benzene','water']
    test_f = Fragmentation(molecule_names)
    
    test_f.molecules = [ SimpleAtoms(name,fragment_list=fragment_names) for name in molecule_names]
    test_f.fragments = [ SimpleAtoms(name) for name in set(fragment_names)]

    test_f.fill_data_reference(data_type='s22')
    test_f.run(write=True)

    error =  0.8573 - test_f.mae
    print 'Test error: {0} eV'.format(error)

#Print "Test load energies"
#test_e = Atomization()
#test_e.add_systems(names)

#test_e.molecules = [SimpleH2Atoms()]
#test_e.atoms = [SimpleH2Atoms()]
#test_e.load_energies('energies.pkl')
#test_e.fill_data_simulated()
#test_e.calculate_error()
#test_e.calculate_mae()
 
        
if __name__ == "__main__":
    main()
