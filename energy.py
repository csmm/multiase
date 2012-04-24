import pickle
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
 #       from ase.data.molecules import molecule, atoms as g2_atoms, g1
        self.molecules = None #List of molecules
        self.atoms = None #Atom object ase-like
        self.calc = None  #Calculator gpaw-like
        self.energies = {}
        self.data = {}
        self.mae = 0.0
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
class SimpleH2Atoms():        
    def __init__(self):
        self.name = 'H2'
        return
    def get_potential_energy(self):
        return 1.0
    def get_chemical_symbols(self):
        return ['H','H']
##atomization GPAW calculator
class GPAWAtoms(name):
    def __init__(self):
        self.name = name
        return
    def set_atoms():
        return
    def set_calculator():
        return
            
names = ['H2']

test = Atomization()
#setting first by hand or calling g22 or s22
test.add_systems(names)
test.molecules = [SimpleH2Atoms()]
test.atoms = [SimpleH2Atoms()]
#calculate
test.get_all_energies()
print test.energies
test.write_energies('energies.pkl')
test.fill_data_simulated()
print test.data
test.calculate_error()
print test.data
test.calculate_mae()
test.write_data('data.pkl')

#to load
#print "test load data"
#test_d = ScalarTest()
#test_d.load_data('data.pkl')
#test_d.calculate_mae()

#print "test load energies"
#test_e = Atomization()
#test_e.add_systems(names)
#test_e.molecules = [SimpleH2Atoms()]
#test_e.atoms = [SimpleH2Atoms()]
#test_e.load_energies('energies.pkl')
#test_e.fill_data_simulated()
#test_e.calculate_error()
#test_e.calculate_mae()

        
    
