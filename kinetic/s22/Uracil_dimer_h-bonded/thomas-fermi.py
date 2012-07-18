#THIS SCRIPT INTEGRATES THE DENSITY TO THE POWER SO THAT THE SCRIPT energy_TF_molecules.py CAN GIVE THE VALUE OF 
#THE KINETIC ENERGY IN THE DIFFERENT APPROXIMATIONS: THOAS-FERMI, GAUSS, BESSEL.

from gpaw import *
from ase import *
from ase.units import *
from numpy import *
from ase.io import write
from ase.io import read
import pickle
from ase.structure import molecule as read_molecule

gridrefinement = 2
Cf = (3.*((3.*(pi**2.))**(2./3.))/10.)
Cg = ((3.*pi)/(2.**(5./3.)))


#Introduce the element of the molecule:
biX = ['Uracil_dimer_h-bonded_f1','Uracil_dimer_h-bonded_f2','Uracil_dimer_h-bonded']

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

def run_per_element(element,read_gpw=False):
     print 'Processing element %s' % element
     fileX = 'data_%s_TF.pkl' % element

#    molecule = read('%s.xyz' % element)
     if read_gpw:
          molecule, calc = restart('%s.gpw' % element)
     else:
          molecule = read_molecule('%s' % element)
          set_work_cell(molecule,h,vacuum)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MOLECULE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#    calc = GPAW(h=0.18, nbands=-4, xc='PBE', occupations=FermiDirac(0.0),fixmom=True)

          calc = GPAW(h=0.18, nbands=-8, xc='PBE')
          if len(molecule)==1: #atom case
               calc.set(hund=True)

          calc.set(occupations=FermiDirac(0.0,fixmagmom=True))
          calc.set(txt='{0}.out'.format(element))
          molecule.set_calculator(calc)
     
          biX_pot = molecule.get_potential_energy()
          calc.write('{0}.gpw'.format(element))

#_______________________________ MOLECULE DENSITY ______________________
#calculate using spin keyword, works for both spin pol and spin unpol

     nmolecule_s0 = 2.*calc.get_all_electron_density(gridrefinement=gridrefinement,pad=False, spin=0)
     nmolecule_s1 = 2.*calc.get_all_electron_density(gridrefinement=gridrefinement,pad=False, spin=1)
     nmolecule_s0 = maximum(nmolecule_s0, 1e-40)  # cuts to avoid numerical errors
     nmolecule_s1 = maximum(nmolecule_s1, 1e-40)  # cuts to avoid numerical errors      
     nmoleculepower_s0 = nmolecule_s0**(5./3.)
     nmoleculepower_s1 = nmolecule_s1**(5./3.)


#_______________________________ MOLECULE INTEGRAL _____________________

     nmoleculebohr_s0 = nmolecule_s0*(Bohr**3)
     nmoleculebohr_s1 = nmolecule_s1*(Bohr**3)

     nmoleculebohrpower_s0 = nmoleculebohr_s0**(5./3.)
     nmoleculebohrpower_s1 = nmoleculebohr_s1**(5./3.)

     
     Imolecule = calc.density.finegd.integrate(nmoleculebohrpower_s0)/2.
     Imolecule += calc.density.finegd.integrate(nmoleculebohrpower_s1)/2.
     KIN_TF = Imolecule*Cf*27.211
     KIN_G = Imolecule*Cg*27.211

     print 'Kinetic energy, TF approx.', KIN_TF
     print 'Kinetic energy, GAUSS approx.', KIN_G

#_______________________________ RESULTS _______________________________

     data1 = {'Molecule density':((nmolecule_s0+nmolecule_s1)/2.),'Molecule density power integral':(Imolecule)}

     output = open(fileX, 'wb')
     pickle.dump(data1, output)
     output.close()
     
for e in biX:
     run_per_element(e,read_gpw=False)
