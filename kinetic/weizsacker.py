from gpaw import *
from ase import *
from ase.units import *
from numpy import *
from ase.io import write
from ase.io import read
import pickle
from ase.structure import molecule as read_molecule



from gpaw.utilities.blas import axpy
from gpaw.fd_operators import Gradient

def calculate_sigma2(calc,n_sg,gd):
    grad_v = [Gradient(gd, v, n=3).apply for v in range(3)]
    gradn_g = gd.empty()
    taut_g = gd.zeros()
    for v in range(3):
        grad_v[v](n_sg, gradn_g)
        axpy(1.0, gradn_g**2, taut_g)
    return taut_g, 1.0

def calculate_sigma(calc,n_sg,gd):
    grad_v = [Gradient(gd, v, n=3).apply for v in range(3)]
    nspins = len(n_sg)
    gradn_svg = gd.empty((nspins, 3))
    sigma_xg = gd.zeros(nspins * 2 - 1)
    for v in range(3):
        for s in range(nspins):
            grad_v[v](n_sg[s], gradn_svg[s, v])
            print'spin, nspins', s, nspins
            axpy(1.0, gradn_svg[s, v] ** 2, sigma_xg[2 * s])
            if nspins == 2:
                axpy(1.0, gradn_svg[0, v] * gradn_svg[1, v],sigma_xg[1])
    return sigma_xg, gradn_svg


gridrefinement = 2

#Introduce the element of the molecule:
#biX = ['O2','HF','H2O','CH4','NH3','CN','CO','F2','HCN','NO']
#biX = ['CO']
biX = ['H']

#Introduce the name of the file for the molecule:
#fileX_w = 'data_%s_W.pkl'

#Introduce:
h = 0.18
vacuum = 4.5

def get_alpha(n, sigma):
        # von Weisaecker
    ind = (n != 0.).nonzero()
    gdms = maximum(sigma, 1e-40)  # |nabla rho|^2
    tau_w = zeros((shape(n)))
    tau_w[ind] = maximum(divide(gdms[ind], 8.0 * n[ind]), 1e-40)
    return tau_w

def set_work_cell(molecule,h,vacuum):
    molecule.center(vacuum=vacuum)
    cell = molecule.get_cell()
    for i in [0,1,2]:
        cell[i,i] = int(cell[i,i]/4./h)*4.*h
    molecule.set_cell(cell)
    return
def run_per_element(element,read_gpw=False):
    print 'Processing element %s' % element
    fileX_w = 'data_%s_W.pkl' % element
    if read_gpw:
        molecule,calc = restart('%s.gpw' %element)
    else:
        try:
            molecule = read_molecule('%s' % element)
        except:
            molecule = read('%s.xyz' % element)
        set_work_cell(molecule,h,vacuum)

        calc = GPAW(h=0.18, nbands=-8, xc='PBE')
        if len(molecule)==1: #atom case
            calc.set(hund=True)

        calc.set(occupations=FermiDirac(0.0,fixmagmom=True))
        calc.set(txt='{0}.out'.format(element))
        molecule.set_calculator(calc)

        biX_pot = molecule.get_potential_energy()
        calc.write('{0}.gpw'.format(element))


    nmolecule_s0 = calc.get_all_electron_density(gridrefinement=gridrefinement,
            pad=False, spin=0)
    nmolecule_s1 = calc.get_all_electron_density(gridrefinement=gridrefinement,
            pad=False, spin=1)

    nmoleculebohr_s0 = 2.*nmolecule_s0*(Bohr**3)
    nmoleculebohr_s1 = 2.*nmolecule_s1*(Bohr**3)

    sigma_s0g, gradn_s0vg = calculate_sigma2(calc,nmoleculebohr_s0,calc.density.finegd)
    sigma_s1g, gradn_s1vg = calculate_sigma2(calc,nmoleculebohr_s1,calc.density.finegd)


    tau_ws0 = get_alpha(nmoleculebohr_s0,sigma_s0g)
    tau_ws1 = get_alpha(nmoleculebohr_s1,sigma_s1g)

    Itau_w = calc.density.finegd.integrate(tau_ws0)/2.
    Itau_w += calc.density.finegd.integrate(tau_ws1)/2.


    print '%%%%%%%%%%%%%%%%%%%% %s %%%%%%%%%%%%%%%%%%' % element
    print 'WEIZSACKER (Ha)', Itau_w
    print 'WEIZSACLER (eV)', Itau_w*27.211
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

    data_w = {'Weizsacker kinetic energy':(Itau_w)}

    output = open(fileX_w, 'wb')
    pickle.dump(data_w, output)
    output.close

for e in biX:
    run_per_element(e, read_gpw=True)

    


