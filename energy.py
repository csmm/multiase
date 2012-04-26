import pickle
import g2_1 as g22_exp_data

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
        self.data_ref = None
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

    def fill_data_reference(self, data_type='g22'):
        for system in self.molecules:
            if data_type == 'g22':
                factor_to_eV=0.0433641146392
                ae = self.ref_data.get_atomization_energy(system.name) * factor_to_eV
            elif data_type == 's22':
                ae = 1.0
            else:
                print "Data type not implemented"
            self.add_data(system.name,ae,'ref')
 
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

    def set_fragments_s22(self,atoms_type='ase'):
        if atoms_type == 'ase':
            f = 0
        elif atoms_type =='test':
            f = [ SimpleAtoms('Benzene',['C','H']), SimpleAtoms('Water',['H','O']) ]
        else:
            print "Atoms type not implemented"
        return f

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
    test.ref_data = g22_exp_data

    test.add_systems(molecule_names)
    test.atoms = [ SimpleAtoms(name,[name]) for name in atom_names]
    test.molecules = [ SimpleAtoms(name,atom_names) for name in molecule_names]

    test.fill_data_reference(data_type='g22')
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

    print "Test fragmentation with s22 set"
    test_f = Atomization()
    molecule_names = ['Benzene-water_complex']
    test_f.add_systems(molecule_names)
    test_f.molecules = [ SimpleAtoms(name,atom_names) for name in molecule_names]
    test_f.atoms = test_f.set_fragments_s22(atoms_type='test')

    test_f.fill_data_reference(data_type='s22')
    test_f.get_all_energies()
    test_f.write_energies('energies.pkl')
    test_f.fill_data_simulated()
    test_f.calculate_error()
    test_f.calculate_mae()
    test_f.write_data('data.pkl')

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
