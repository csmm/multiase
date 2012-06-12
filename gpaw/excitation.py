from ase import Atoms, Atom
from ase.io import read, write
from ase.units import Bohr
from gpaw import GPAW, restart
from gpaw.lrtddft import LrTDDFT
from gpaw import dscf
import os

def write_gpaw(atoms,calc,outfile):
    atoms.set_calculator(calc)
    e = atoms.get_potential_energy()
    calc.write(outfile,mode='all')
    return 
    
def optical(gpwfile,istart,jend):
    """istart= band index of the first occ. band to consider
    jend band index of the last unocc. band to consider"""
    atoms, calc = restart(gpwfile)
    lr = LrTDDFT(calc, xc='LDA', istart=istart, jend=jend,
                 nspins=1) # only singlet
    lr.write('lr.dat.gz')
    return

def write_bader(atoms,calc,outfile):
    n = calc.get_all_electron_density(gridrefinement=2)
    n *= Bohr**3
    write(outfile, atoms, data=n)
    return
def bader_charge(gpwfile,outfile):
    atoms, calc = restart(gpwfile)
    write_bader(atoms,calc,outfile)
    return
def write_deltascf(outfile,n_homo,nbands,h):
    atoms, calc = restart(outfile)
    calc.set(nbands=nbands, h=h, xc='PBE', spinpol=False,
             convergence={'energy': 100,
                          'density': 100,
                          'eigenstates': 1.0e-9,
                          'bands': nbands+1})

    E_gs = atoms.get_potential_energy()

    # Obtain the pseudowavefunctions and projector overlaps of the
    n = n_homo+1
    molecule = range(len(atoms))
    wf_u = [kpt.psit_nG[n] for kpt in calc.wfs.kpt_u]
    p_uai = [dict([(molecule[a], P_ni[n]) for a, P_ni in kpt.P_ani.items()])
             for kpt in calc.wfs.kpt_u]

    n = n_homo
    wf_u_homo = [kpt.psit_nG[n] for kpt in calc.wfs.kpt_u]
    p_uai_homo = [dict([(molecule[a], P_ni[n]) for a, P_ni 
                        in kpt.P_ani.items()])
                  for kpt in calc.wfs.kpt_u]

    # Excited state calculation
    #--------------------------------------------
    atoms_es = atoms.copy()
    calc_es = GPAW(nbands=nbands, h=h, xc='PBE', spinpol=False,
                   convergence={
            'energy': 100,
            'density': 100,
            'eigenstates': 1.0e-9,
            'bands': 'occupied'},
                   maxiter=1000)

    atoms_es.set_calculator(calc_es)
    lumo = dscf.AEOrbital(calc_es, wf_u, p_uai)
    homo =  dscf.AEOrbital(calc_es, wf_u_homo, p_uai_homo)
    dscf.dscf_calculation(calc_es, [[1.0,lumo,0]]
                          , atoms_es)
    
    E_es = atoms_es.get_potential_energy()
    
    from gpaw.utilities.dscftools import dscf_collapse_orbitals 
    dscf_collapse_orbitals(calc_es,verify_density=True)
    calc_es.print_all_information()

    calc_es.write('dscf.gpw',mode='all')
    return


## atoms
atoms = read('relax-mm-PBE.traj')
atoms.set_initial_magnetic_moments([0.0,0.0])
## gpaw calculation
h=0.3
nbands=-2
charge = 0
calc = GPAW(h=h, nbands=nbands, xc='PBE', txt='out.txt',charge=charge,
            convergence={'bands':nbands+1}, spinpol=False)

gpwfile = 'out.gpw'
write_gpaw(atoms,calc,gpwfile)

n_homo = 10
istart = n_homo-2
jend = n_homo+2
optical(gpwfile,istart,jend)

bader_charge(gpwfile,'density.cube')

write_deltascf(gpwfile,n_homo, nbands,h)
bader_charge('dscf.gpw','density_dsc.cube')
