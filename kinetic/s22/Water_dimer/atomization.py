#THIS SCRIPT CALCULATES THE ATOMIZATION ENERGY OF MOLECULES (organic)

from gpaw import *
from ase import *
from ase.units import *
from numpy import *
from ase.io import write
from ase.io import read
import pickle
from ase.structure import molecule as read_molecule
import numpy as np

gridrefinement = 2
Cf = (3.*((3.*(pi**2.))**(2./3.))/10.)
Cg = ((3.*pi)/(2.**(5./3.)))

sumW = 0.
sumTF = 0.
sumG = 0.
sumC = 0.
sumCnew = 0.

data_sum = {'sumW':(sumW),'sumTF':(sumTF),'sumG':(sumG),'sumC':(sumC),'sumCnew':(sumCnew)}
sums = np.zeros(5)

#Introduce:
h = 0.18
vacuum = 4.5

#Introduce the molecules to be studied:
MOL = ['Water_dimer']

#Introduce the different atoms that are in the molecules:
ATO = ['H','H','O','H','H','O']

#%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUM OF ATOM ENERGIES %%%%%%%%%%%%%%%%%%%%%%%%%

#output = open('sums.pkl', 'wb')
#pickle.dump(data_sum, output)
#output.close()


def run_per_atom(sums, element, sumW, sumTF, sumG, sumC, sumCnew):

    print 'Processing element %s' % element

    fileX = 'data_%s_TF.pkl' % element
    fileX_w = 'data_%s_W.pkl' % element
    fileX_C = 'data_%s_C.pkl' % element
    fileX_Cnew = 'data_%s_Cnew.pkl' % element

    #Getting Weizsacker kinetic energy.
    pkl_file = open(fileX_w, 'rb')
    data_w = pickle.load(pkl_file)
    pkl_file.close

    Itau_w = data_w['Weizsacker kinetic energy']
    Wkin = Itau_w

    #Getting Thomas-Fermi kinetic energy.
    pkl_file = open(fileX, 'rb')
    data1 = pickle.load(pkl_file)
    pkl_file.close()

    Iatom = data1['Molecule density power integral']
    TFkin = Iatom*Cf

    #Getting Gauss kinetic energy.
    Gkin = Iatom*Cg

    #Getting combined kinetic energy.
    pkl_file = open(fileX_C, 'rb')
    data_c = pickle.load(pkl_file)
    pkl_file.close()

    Ckin = data_c['Combined kinetic energy']

    #Getting new combined kinetic energy.
    pkl_file = open(fileX_Cnew, 'rb')
    data_cnew = pickle.load(pkl_file)
    pkl_file.close

    Cnewkin = data_cnew['New combined kinetic energy']

    sums[0] += Wkin
    sums[1] += TFkin
    sums[2] += Gkin
    sums[3] += Ckin
    sums[4] += Cnewkin

for e in ATO:
    run_per_atom(sums, e, sumW, sumTF, sumG, sumC, sumCnew)
    
    data_sum = {'sumW':(sums[0]),'sumTF':(sums[1]),'sumG':(sums[2]),'sumC':(sums[3]),'sumCnew':(sums[4])}
    
    output = open('sums.pkl', 'wb')
    pickle.dump(data_sum, output)
    output.close()


#%%%%%%%%%%%%%%%%%%%%%%%%%% MOLECULE ENERGY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Wkin_mol = 0.
TFkin_mol = 0.
Gkin_mol = 0.
Ckin_mol = 0.
Cnewkin_mol = 0.

def run_per_molecule(element, Wkin_mol, TFkin_mol, Gkin_mol, Ckin_mol, Cnewkin_mol):
    print 'Processing molecule %s' % element

    fileX = 'data_%s_TF.pkl' % element
    fileX_w = 'data_%s_W.pkl' % element
    fileX_C = 'data_%s_C.pkl' % element
    fileX_Cnew = 'data_%s_Cnew.pkl' % element

    #Getting Weizsacker kinetic energy.
    pkl_file = open(fileX_w, 'rb')
    data_w = pickle.load(pkl_file)
    pkl_file.close

    Itau_w = data_w['Weizsacker kinetic energy']
    Wkin_mol = Itau_w

    #Getting Thomas-Fermi kinetic energy.
    pkl_file = open(fileX, 'rb')
    data1 = pickle.load(pkl_file)
    pkl_file.close()

    Iatom = data1['Molecule density power integral']
    TFkin_mol = Iatom*Cf

    #Getting Gauss kinetic energy.
    Gkin_mol = Iatom*Cg

    #Getting combined kinetic energy.
    pkl_file = open(fileX_C, 'rb')
    data_c = pickle.load(pkl_file)
    pkl_file.close()

    Ckin_mol = data_c['Combined kinetic energy']

    #Getting new combined kinetic energy.
    pkl_file = open(fileX_Cnew, 'rb')
    data_cnew = pickle.load(pkl_file)
    pkl_file.close()

    Cnewkin_mol = data_cnew['New combined kinetic energy']

    data_mol = {'Wkin_mol':(Wkin_mol),'TFkin_mol':(TFkin_mol),'Gkin_mol':(Gkin_mol),'Ckin_mol':(Ckin_mol),'Cnewkin_mol':(Cnewkin_mol)}

    output = open('mol_%s.pkl' % element, 'wb')
    pickle.dump(data_mol, output)
    output.close()


for i in MOL:
    run_per_molecule(i, Wkin_mol, TFkin_mol, Gkin_mol, Ckin_mol, Cnewkin_mol)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ATOMIZATION ENERGY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def run_per_element(element):

    pkl_file = open('sums.pkl', 'rb')
    data_sum = pickle.load(pkl_file)
    pkl_file.close()

    sumTF = data_sum['sumTF']
    sumG = data_sum['sumG']
    sumW = data_sum['sumW']
    sumC = data_sum['sumC']
    sumCnew = data_sum['sumCnew']

    pkl_file = open('mol_%s.pkl' % element, 'rb')
    data_mol = pickle.load(pkl_file)
    pkl_file.close()

    TFkin_mol = data_mol['TFkin_mol']
    Gkin_mol = data_mol['Gkin_mol']
    Wkin_mol = data_mol['Wkin_mol']
    Ckin_mol = data_mol['Ckin_mol']
    Cnewkin_mol = data_mol['Cnewkin_mol']

    print ' ' 
    print 'sumTf', sumTF*27.211
    print 'sumG', sumG*27.211
    print 'sumW', sumW*27.211
    print 'sumC', sumC*27.211
    print 'sumCnew', sumCnew*27.211
    print ' '
    print 'TFkin_mol', TFkin_mol*27.211
    print 'Gkin_mol', Gkin_mol*27.211
    print 'Wkin_mol', Wkin_mol*27.211
    print 'Ckin_mol', Ckin_mol*27.211  
    print 'Cnewkin_mol', Cnewkin_mol*27.211
    
    
    print ' '
    print '%%%%%%%%%%%% ATOMIZATION ENERGY OF %s MOLECULE %%%%%%%%%%%%' % element
    print ' '
    print 'THOMAS-FERMI approx:  ', (sumTF - TFkin_mol)*27.211, 'eV'
    print 'GAUSS approx:         ', (sumG - Gkin_mol)*27.211, 'eV'
    print 'WEIZSACKER approx:    ', (sumW - Wkin_mol)*27.211, 'eV'
    print 'COMBINED approx:      ', (sumC - Ckin_mol)*27.211, 'eV'
    print 'NEW COMBINED approx:  ', (sumCnew - Cnewkin_mol)*27.211, 'eV'
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

for a in MOL:
    run_per_element(a)















































