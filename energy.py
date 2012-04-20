import pickle
class ScalarTest(object):
    def __init__(self):
        self.data = {}
        self.mae = 0.0
    def add_systems(self, list_names): #it will be expanded to dyn error
        for name in list_names:
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
        print 'Min absolute error is {}'.format(self.mae)
    def write_data(self,filename):
        f = open(filename,'wb')
        pickle.dump(self.data, f)
        f.close()
        f = open('data.txt','w')
        for k,v in self.data.iteritems():
            f.write('{} {}\n'.format(k,v))
        f.write('mae: {}'.format(self.mae))
    def load_data(self,filename):
        f = open(filename,'rb')
        self.data = pickle.load(f)
        f.close()

class Atomization(ScalarTest):
    def __init__(self):
 #       from ase.data.molecules import molecule, atoms as g2_atoms, g1

        self.molecules = None
        self.atoms = None
        self.calc = None
        self.energies = {}

    def get_all_energies(self):
        for atom in self.atoms:
            e = atom.get_potential_energy()
            self.energies[atom.get_chemical_symbols[0]] = e
        for system in self.molecules:
            e = system.get_potential_energy()
            self.energies[system.name] = e
        def fill_data(self):
            for system in self.molecules:
                ae = -1.0*energies[system.name]
                for atom in system.get_chemical_symbols():
                    ae += energies[atom]
                add_data(system,ae,'sim')
            
names = ['H2']

test = Atomization()
test.add_systems(names)
test.calculate_error()
test.write_data('data.pkl')

#to load
#test_new = AccuracyTest()
#test_new.load_data('data.pkl')
#test_new.calculate_mae()

        
    
