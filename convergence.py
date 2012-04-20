from ase import Atoms, Atom
from ase.io import read
from gpaw import GPAW
from ase.parallel import parprint
from math import fabs

def atomization_energy(atom_list,molecule_name,h,cell):

# molecule:
    molecule = read(molecule_name+'.xyz')
    molecule.set_cell(cell)
    molecule.center()

    calc = GPAW(h=h, nbands=-5, xc='PBE', txt=molecule_name+'.out',spinpol=True)
    molecule.set_calculator(calc)

    e2 = molecule.get_potential_energy()
    #calc.write(molecule_name+'.gpw')

# atoms:
    e1 = 0.0
    for atom_name,n_a in atom_list:
        atom = Atoms(atom_name,
                     positions=[(0, 0, 0)],
                     cell=cell)
        atom.center()

        atom.set_calculator(calc)
        calc.set(txt=atom_name+'.out',hund=True)
        e1 += n_a*atom.get_potential_energy()
        #calc.write(atom_name+'.gpw')


    parprint('\t\tatom energy:     %5.2f eV' % e1)
    parprint('\t\tmolecule energy: %5.2f eV' % e2)
    parprint('\t\tatomization energy:       %5.2f eV' % (e1 - e2))
    return e1-e2

def get_work_cell(molecule_file,h,vacuum):
    molecule = read(molecule_file)
    molecule.center(vacuum=vacuum)
    cell = molecule.get_cell()
    for i in [0,1,2]:
        cell[i,i] = int(cell[i,i]/4./h)*4.*h
    return cell

def run_vacuum(vacuum,h,molecule_name):
    parprint('vacuum:',vacuum)
    cell = get_work_cell(molecule_name+'.xyz',h,vacuum)
    return atomization_energy(atom_list,molecule_name,h,cell)

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

#######PARAMETERS##############
atom_list = [('H',2)] # [ ('atom1',number),('atom2',number),...]
molecule_name = 'H2'
tol = 1.00 #eV

v_init = 3.0 #Ang
dv = 0.5 
v_max = 6.0

h_init = 0.4 #
dh = -0.1
h_min = 0.2
#########################

h = h_init
conv = False
ae = 1000

while conv == False:
    parprint('h:',h)
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



