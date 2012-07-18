#This script will provide with the kinetic energy calculated in the combined approximation. The combination is between Weizsacker and Thomas Fermi. 

from gpaw import *
from ase import *
from ase.units import *
from numpy import *
import pickle

#Thomas-Fermi constant.
Cf = (3.*((3.*(pi**2.))**(2./3.))/10.)

#Introduce the molecules to be studied:
biX = ['Adenine-thymine_Watson-Crick_complex_f1','Adenine-thymine_Watson-Crick_complex_f2','Adenine-thymine_Watson-Crick_complex']

#Introduce the number of electrons in the atom or molecule:
N = [70.,66.,136.]

def run_per_element(element, N):
    print ' '
    print 'processing element %s' % element

    fileX = 'data_%s_TF.pkl' % element
    fileX_w = 'data_%s_W.pkl' % element
    fileX_C = 'data_%s_C.pkl' % element
    
    #Getting Weizsacker kinetic energy.
    pkl_file = open(fileX_w, 'rb')
    data_w = pickle.load(pkl_file)
    pkl_file.close

    Itau_w = data_w['Weizsacker kinetic energy']
    Wkin = Itau_w

    #Getting Thomas-Fermi kinetic energy.
    pkl_file = open(fileX,'rb')
    data1 = pickle.load(pkl_file)
    pkl_file.close()

    Iatom = data1['Molecule density power integral']
    TFkin = Iatom*Cf


    #Getting the factor that depends on N.
    #N = read('%s' % element)

    gamma = 1.0 - (1.412/(N**(1./3.))) 
    
    com_kin = Wkin + (gamma*TFkin)

    print ' '
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print 'Gamma:', gamma
    print ' '
    print 'COMBINED KINETIC ENERGY, %s:' % element, com_kin*27.211
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print ' '

    data_c = {'Combined kinetic energy':(com_kin)}

    output = open(fileX_C, 'wb')
    pickle.dump(data_c, output)
    output.close

for e,n in zip(biX,N):
    run_per_element(e,n)

