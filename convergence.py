from ase import Atoms, Atom
from ase.io import read
from ase.structure import molecule as mol
from gpaw import GPAW
from ase.parallel import parprint
from math import fabs

def atomization_energy(molecule_name,h,cell,k=1,pbc=False):

# molecule:
#    try:
    print molecule_name
    molecule = mol(molecule_name)
#    except:
#        molecule = read(molecule_name+'.xyz')
    molecule.set_cell(cell)
    if pbc:
        molecule.set_pbc(True)
    else:
        molecule.center()
    print (k,k,k)    
    calc = GPAW(h=h, nbands=-5, xc='PBE', txt=molecule_name+'.out',
                spinpol=True, kpts=(k, k, k), usesymm=False)
    molecule.set_calculator(calc)

    e2 = molecule.get_potential_energy()
    #calc.write(molecule_name+'.gpw')

# atoms:
    e1 = 0.0
    atom_list = molecule.get_chemical_symbols()
    for atom_name in set(atom_list):
        atom = Atoms(atom_name,
                     positions=[(0, 0, 0)],
                     cell=cell,
                     pbc=pbc )
        atom.center()

        atom.set_calculator(calc)
        calc.set(txt=atom_name+'.out',hund=True)
        n_a = atom_list.count(atom_name)
        e1 += n_a*atom.get_potential_energy()
        #calc.write(atom_name+'.gpw')


    parprint('\t\tatoms energy:     %5.2f eV' % e1)
    parprint('\t\tmolecule energy: %5.2f eV' % e2)
    parprint('\t\tatomization energy:       %5.2f eV' % (e1 - e2))
    return e1-e2

def get_work_cell(molecule_name,h,vacuum):
#    molecule = read(molecule_file+'.xyz')
    molecule = mol(molecule_name)
    molecule.center(vacuum=vacuum)
    cell = molecule.get_cell()
    for i in [0,1,2]:
        cell[i,i] = int(cell[i,i]/4./h)*4.*h
    return cell

def run_vacuum(vacuum,h,molecule_name):
    parprint('vacuum:',vacuum)
    cell = get_work_cell(molecule_name,h,vacuum)
    return atomization_energy(molecule_name,h,cell)

def run_k(k,h,cell,molecule_name):
    parprint('k:',k)
    #h = get_work_cell(molecule_name+'.xyz',h,vacuum)
    return atomization_energy(molecule_name,h,cell,pbc=True,k=k)


def converge_vacuum(h,v_init,dv,vmax,tol):
    conv = False
    vacuum = v_init
    ae = 10000
    while conv == False:
        ae_new = run_vacuum(vacuum,h,molecule_name)
        if fabs(ae-ae_new) >= tol:
            vacuum += dv
            ae = ae_new
            if vacuum > vmax:
                parprint('vacuum not converged')
                stop
        else:
            parprint('box size converged, sugggested vacuum (ang)',vacuum-dv)
            conv = True
    return ae


def converge_kpoints(h,cell,k_init,dk,kmax,tol):
    conv = False
    k = k_init
    ae = 10000
    while conv == False:
        ae_new = run_k(k,h,cell,molecule_name)
        if fabs(ae-ae_new) >= tol:
            k += dk
            ae = ae_new
            if k > kmax:
                parprint('kpoints not converged')
                stop
        else:
            parprint('kpoints converged, sugggested k',k-dk)
            conv = True
    return ae

#######PARAMETERS##############
#atom_list = [('H',2)] # [ ('atom1',number),('atom2',number),...]
molecule_name = 'H2'
tol = 10.0 #eV

v_init = 3.5 #Ang
dv = 0.5 
v_max = 6.0

h_init = 0.3 #
dh = -0.1
h_min = 0.16

k_test = False
k_init = 1
k_max = 20
dk = 1
cell = [4.0,4.0,4.0]
#########################

h = h_init
conv = False
ae = 1000

while conv == False:
    parprint('h:',h)
    if k_test:
        ae_new= converge_kpoints(h,cell,k_init,dk,k_max,tol)
    else:
        ae_new = converge_vacuum(h,v_init,dv,v_max,tol)
        
    if fabs(ae-ae_new) >= tol:
        h += dh
        ae = ae_new
        if h < h_min:
            parprint('grid not converged')
            stop
    else:
        parprint('grid converged, suggested grid (Ang)', h-dh)
        parprint('tolerance used (eV)', tol)
        conv = True



