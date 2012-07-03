from gpaw import *
from ase import *
from ase.units import *
from numpy import *
from ase.io import write
from ase.io import read
import pickle
from ase.structure import molecule as read_molecule


from gpaw.xc.lda import LDA
from gpaw.utilities.blas import axpy
from gpaw.fd_operators import Gradient
#from gpaw.sphere.lebedev import Y_nL, weight_n
#from gpaw.xc.pawcorrection import rnablaY_nLv

biX = ['CH4','NH3','H2O','HF','CN','HCN','CO','N2','O2','F2']



#def set_grid_descriptor(self, gd):
    #LDA.set_grid_descriptor(self, gd)

def calculate_sigma2(calc,n_sg,gd):
    grad_v = [Gradient(gd, v, n=3).apply for v in range(3)]
    #nspins = calc.wfs.nspins
    nspins = 1
    #taut_sG = calc.wfs.gd.empty((nspins))
    gradn_g = gd.empty()
    print 'gradn_g', shape(gradn_g)
    print 'n_sg', shape(n_sg)
    #for s in range(nspins):
    taut_g = gd.zeros()
    for v in range(3):
        grad_v[v](n_sg, gradn_g)
        axpy(1.0, gradn_g**2, taut_g)
      #  ntinv_g = 0. * taut_g
      #  nt_ok = np.where(n_sg[s] > 1e-7)
      #  ntinv_g[nt_ok] = 1.0 / nt_sg[s][nt_ok]
      #  taut_g *= ntinv_g
      #  self.restrict(taut_g, taut_sG[s])
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


gridrefinement = 1

#Introduce the element of the molecule:
#biX = ['O2','HF','H2O','CH4','NH3','CN','CO','F2','HCN','NO']

#biX = ['O']

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
         molecule.set_cell(cell)def run_per_element(element):
    print 'Processing element %s' % element

    fileX_w = 'data_%s_W.pkl' % element

    try:
        molecule = read_molecule('%s' % element)
    except:
        molecule = read('%s.xyz' % element)
    set_work_cell(molecule,h,vacuum)


    calc = GPAW(h=0.18, nbands=-8, xc='PBE',
            occupations=FermiDirac(0.0,fixmagmom=True))


    calc.set(txt='{0}.out'.format(element))
    molecule.set_calculator(calc)

    biX_pot = molecule.get_potential_energy()
    calc.write('{0}.gpw'.format(element))


    nmolecule = calc.get_all_electron_density(gridrefinement=gridrefinement,
            pad=False)
    nmoleculebohr = nmolecule*(Bohr**3)

    sigma_xg, gradn_svg = calculate_sigma2(calc,nmoleculebohr,calc.density.gd)


    tau_w = get_alpha(nmoleculebohr,sigma_xg)

    Itau_w = calc.density.gd.integrate(tau_w)

    print '%%%%%%%%%%%%%%%%%%%% %s %%%%%%%%%%%%%%%%%%' % element
    print 'WEIZSACKER (Ha)', Itau_w
    print 'WEIZSACLER (eV)', Itau_w*27.211
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

    data_w = {'Weizsacker kinetic energy':(Itau_w)}

    output = open(fileX_w, 'wb')
    pickle.dump(data_w, output)
    output.close

for e in biX:
    run_per_element(e)

         return


