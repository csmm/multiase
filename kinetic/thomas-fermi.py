#THIS SCRIPT INTEGRATES THE DENSITY TO THE POWER SO THAT THE SCRIPT energy_TF_molecules.py CAN GIVE THE VALUE OF THE KINETIC ENERGY IN THE DIFFERENT APPROXIMATIONS: THOAS-FERMI, GAUSS, BESSEL.

from gpaw import *
from ase import *
from ase.units import *
from numpy import *
from ase.io import write
from ase.io import read
import pickle
from ase.structure import molecule as read_molecule

gridrefinement = 1

#Introduce the element of the molecule:
biX = ['O2','HF','H2O','CH4','NH3','CN','CO','F2','HCN','N2']

#Introduce the name of the file for the molecule:
#fileX = 'data_%s_TF.pkl'

#Introduce:
h = 0.18
vacuum = 4.5

def set_work_cell(molecule,h,vacuum):
     molecule.center(vacuum=vacuum)
     cell = molecule.get_cell()
     for i in [0,1,2]:
        cell[i,i] = int(cell[i,i]/4./h)*4.*h
     molecule.set_cell(cell)
return

def run_per_element(element):
    print 'Processing element %s' % element
#    molecule = read('%s.xyz' % element)
    molecule = read_molecule('%s' % element)
    set_work_cell(molecule,h,vacuum)

    fileX = 'data_%s_TF.pkl' % element



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MOLECULE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#    calc = GPAW(h=0.18, nbands=-4, xc='PBE', occupations=FermiDirac(0.0),fixmom=True)

    calc = GPAW(h=0.18, nbands=-4, xc='PBE',
            occupations=FermiDirac(0.0,fixmagmom=True))


    calc.set(txt='{0}.out'.format(element))
    molecule.set_calculator(calc)

    biX_pot = molecule.get_potential_energy()
    calc.write('{0}.gpw'.format(element))

#_______________________________ MOLECULE DENSITY ______________________

    nmolecule = calc.get_all_electron_density(gridrefinement=gridrefinement)

    nmoleculepower = nmolecule**(5./3.)

#_______________________________ MOLECULE INTEGRAL _____________________

    nmoleculebohr = nmolecule*(Bohr**3)

    nmoleculebohrpower = nmoleculebohr**(5./3.)

    Imolecule = calc.density.gd.integrate(nmoleculebohrpower)

#_______________________________ RESULTS _______________________________

    data1 = {'Molecule density':(nmolecule),'Molecule density power integral':(Imolecule)}

    output = open(fileX, 'wb')
    pickle.dump(data1, output)
    output.close()

for e in biX:
    run_per_element(e)
