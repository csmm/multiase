import pickle
from ase.structure import molecule as mol
from gpaw import GPAW
from gpaw.occupations import FermiDirac
#from ase.data.molecules import molecule, atoms as g2_atoms, g1
from ase.data.g2_1 import atom_names as g22_atom_names 
from ase.data.g2_1 import molecule_names as g22_molecule_names 
from ase.data.g2_1 import data as g22_data
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
    def fill_data_experimental(self):   
        for system in self.molecules:
            ae = -1.0 * self.exp_data[system.name]['thermal correction']
            for atom in system.get_chemical_symbols():
                ae += self.exp_data[atom]['thermal correction']
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
class SimpleH2Atoms():        
    def __init__(self):
        self.name = 'H2'
        return
    def get_potential_energy(self):
        return 1.0
    def get_chemical_symbols(self):
        return ['H','H']
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
        else:
            e = 1000 #not converged value
        return e
    def get_chemical_symbols(self):
        return self.system.get_chemical_symbols()


test = Atomization()
#calling g22 or s22
molecule_names = g22_molecule_names # ['H2,...]
atom_names = g22_atom_names #['H',...]
test.exp_data = g22_data

test.add_systems(molecule_names)
test.atoms = [GPAWSystems(name,h=0.17,vacuum=6.0,xc='PBE') 
              for name in atom_names]
test.molecules = [GPAWSystems(name,h=0.4,vacuum=6.0,xc='PBE') 
                  for name in molecule_names]

test.fill_data_experimental()
print test.data
stop
#calculate
test.get_all_energies()
test.write_energies('energies.pkl')
test.fill_data_simulated()
test.calculate_error()
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

        
    
