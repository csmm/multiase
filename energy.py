import pickle
import g2_1 as exp_data

class ScalarTest(object):
    def __init__(self):
        self.data = {}
        self.mae = 0.0
    def add_systems(self,names): #it will be expanded to dyn error
        for name in names:
            self.data[name] = {'exp':0.0,'sim':0.0,'err':0.0}
    def add_data(self, name, data, data_type):
        """ Data types are 'exp' or 'sim' """
        self.data[name][data_type] = data
    def calculate_error(self):
        for name, values in self.data.iteritems():
            values['err'] = values['exp']-values['sim']
    def calculate_mae(self):
        for name, values in self.data.iteritems():
            self.mae += abs(values['exp']-values['sim'])
        self.mae /= len(self.data)
        print 'Mean absolute error is {0} eV'.format(self.mae)
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

class Atomization(ScalarTest):
    def __init__(self):
        self.molecules = None #List of molecules
        self.atoms = None #Atom object ase-like
        self.calc = None  #Calculator gpaw-like
        self.energies = {}
        self.data = {}
        self.mae = 0.0
        self.exp_data = None
    def get_all_energies(self):
        for atom in self.atoms:
            e = atom.get_potential_energy()
            self.energies[atom.get_chemical_symbols()[0]] = e
        for system in self.molecules:
            e = system.get_potential_energy()
            self.energies[system.name] = e
    
    def fill_data_simulated(self):
        for system in self.molecules:
            ae = -1.0 * self.energies[system.name]
            for atom in system.get_chemical_symbols():
                ae += self.energies[atom]
            self.add_data(system.name,ae,'sim')

    def fill_data_experimental(self, factor_to_eV=1.0):
        for system in self.molecules:
            ae = self.exp_data.get_atomization_energy(system.name) * factor_to_eV
            self.add_data(system.name,ae,'exp')

    def write_energies(self,filename):
        f = open(filename,'wb')
        pickle.dump(self.energies, f)
        f.close()
        f = open('energies.txt','w')
        for k,v in self.energies.iteritems():
            out = str(k)+'\t'+str(v)+'\n'
            f.write(out)
#        f.write('mae: {0}'.format(self.mae))
    def load_energies(self,filename):
        f = open(filename,'rb')
        self.energies = pickle.load(f)
        f.close()

##simple atoms for testing
class SimpleAtoms():        
    def __init__(self,name,chemical_symbols):
        self.name = name
        self.chemical_symbols = chemical_symbols
        return
    def get_potential_energy(self):
        return 1.0
    def get_chemical_symbols(self):
        return self.chemical_symbols

def main():
    print "Test atomization"
    test = Atomization()
    molecule_names = ['CH4']
    atom_names = ['H','C']
    test.exp_data = exp_data

    test.add_systems(molecule_names)
    test.atoms = [ SimpleAtoms(name,[name]) for name in atom_names]
    test.molecules = [ SimpleAtoms(name,atom_names) for name in molecule_names]

    test.fill_data_experimental(factor_to_eV=0.0433641146392)
    test.get_all_energies()
    test.write_energies('energies.pkl')
    test.fill_data_simulated()
    test.calculate_error()
    test.calculate_mae()
    test.write_data('data.pkl')

    error = 17.2206599701 - test.mae
    print 'Test error: {0} eV'.format(error)

    #to load
    print "Test load data"
    test_d = ScalarTest()
    test_d.load_data('data.pkl')
    test_d.calculate_mae()
    error = 17.2206599701 - test.mae
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
