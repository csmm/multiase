#This script will provide with the kinetic energy calculated in the combined approximation. The combination is between Weizsacker and Thomas Fermi. 

from gpaw import *
from ase import *
from ase.units import *
from numpy import *
import pickle

#Thomas-Fermi constant.
Cf = (3.*((3.*(pi**2.))**(2./3.))/10.)

#Introduce the molecules to be studiedbiX = ['CH4']
biX = ['F2']

#N_H = 1.
#N_B = 5.
#N_C = 6.
#N_N = 7.
#N_O = 8.
#N_F = 9.

#Introduce the name of the file of the molecule in TF:
#fileX = 'data_%s_TF.pkl' 

#Introduce the name of the file of the molecule in W:
#fileX_w = 'data_%s_W.pkl' 



def run_per_element(element):
    print ' '
    print 'processing element %s' % element

    fileX = 'data_%s_TF.pkl' % element
    fileX_w = 'data_%s_W.pkl' % element

    #Getting Thomas-Fermi kinetic energy.
    #pkl_file = open(fileX,'rb')
    #data1 = pickle.load(pkl_file)
    #pkl_file.close()

    #Iatom = data1['Atom density power integral']
    #TFkin = Iatom*Cf

    #Getting Weizsacker kinetic energy.
    pkl_file = open(fileX_w, 'rb')
    data_w = pickle.load(pkl_file)
    pkl_file.close

    Itau_w = data_w['Weizsacker kinetic energy']
    Wkin = Itau_w*27.211

    #Getting Thomas-Fermi kinetic energy.
    pkl_file = open(fileX,'rb')
    data1 = pickle.load(pkl_file)
    pkl_file.close()

    Iatom = data1['Molecule density power integral']
    TFkin = Iatom*Cf


    #Getting the factor that depends on N.
    #N = read('%s' % element)
    N = 18.

    gamma = 1.0 - (1.412/(N**(1./3.))) 
    
    com_kin = Wkin + (gamma*TFkin)

    print ' '
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print 'COMBINED KINETIC ENERGY, %s:' % element, com_kin
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print ' '


for e in biX:
    run_per_element(e)

