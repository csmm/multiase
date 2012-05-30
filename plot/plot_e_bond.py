#red blue green rbg, dashes squares triangles -- s ^
from ase.io import read
from ase.io.trajectory import PickleTrajectory
from ase.calculators.neighborlist import NeighborList
import matplotlib.pyplot as plt
from ase.units import fs, kB
from matplotlib.backends.backend_pdf import PdfPages

def plot_energies(path,md_type,title,t0_histo,bin_histo,write=False):
    traj = PickleTrajectory(path+md_type+'.traj')
    ref_e = traj[0].get_total_energy()
    ref_epot = traj[0].get_potential_energy()
    e = [ atoms.get_total_energy() - ref_e for atoms in traj]
    ekin = [ atoms.get_kinetic_energy() for atoms in traj] 
    epot = [ atoms.get_potential_energy() - ref_epot for atoms in traj]
    temp = [ atoms.get_temperature() for atoms in traj]

    plt.figure(1)
    plt.title(title)
    plt.ylabel('Energy [eV]')
    plt.xlabel('Time [fs]')
    plt.plot( e, 'g-' , label=r'Etot[t]-Etot[0]')
    plt.legend(loc='upper right')
    if write:
        plt.savefig(md_type+'_etot.eps', format='eps')
        plt.savefig(pp, format='pdf')
    plt.figure(2)
    plt.subplot(211)
    plt.title(title)
    plt.ylabel('Energy [eV]')
    plt.xlabel('Time [fs]')
    plt.plot( epot, 'g-' , label=r'Epot[t]-Epot[0]')
    plt.legend(loc='upper right')

    plt.subplot(212)
    plt.ylabel('Energy [eV]')
    plt.xlabel('Time [fs]')
    plt.plot( ekin, 'b-' ,label=r'Ekin')
    plt.legend(loc='upper right')
    if write:
        plt.savefig(md_type+'_e.eps', format='eps')
        plt.savefig(pp, format='pdf')
    plt.figure(3)
    plt.title(title)
    plt.subplot(211)
    plt.xlabel('Temperature [K]')
    plt.annotate("using t0 = {0} fs".format(t0_histo), xy=(0.75, 0.75),
                 xycoords="axes fraction")
    plt.hist(temp[t0_histo:], bin_histo)
    
    plt.subplot(212)
    plt.ylabel('Temperature [K]')
    plt.xlabel('Time [fs]')
    plt.plot(temp,  'b-' )
    if write:
        plt.savefig(md_type+'_T.eps', format='eps')
        plt.savefig(pp, format='pdf')
    return

def plot_oh(path,md_type,title,t0_histo,bin_histo,write=False):
    traj = PickleTrajectory(path+md_type+'.traj')
    oh = []
    oh_distances = []
    o_index = [atom.index for atom in traj[0] if atom.symbol == 'O']
    h_index = [atom.index for atom in traj[0] if atom.symbol == 'H']
    bonds_init = [ (o,h,traj[0].get_distance(o,h,mic=True))
                   for o in o_index for h in h_index]
    bonds_pva = [ (o,h) for o,h,d in bonds_init if d < 1.1 ]
    no = len(o_index)
    assert len(bonds_pva) == no

    for atoms in traj:
        bonds = [ (o,h,atoms.get_distance(o,h,mic=True)) 
                  for o1,h in bonds_pva for o in o_index]
        bonds_pva_traj = [ d for o,h,d in bonds if d < 1.2 ]
        assert len(bonds_pva_traj) == no
        bonds_oh = [ d for o,h,d in bonds if d <= 2.0]
        bonds_no_pva =  [ d for d in bonds_oh if d > 1.1]

        oh  += [1.0*len(bonds_no_pva)/no]
        oh_distances += bonds_no_pva
        
    plt.figure(4)
    plt.subplot(211)
    plt.title(title)
    plt.ylabel('Number of bonds/O atom')
    plt.xlabel('Time [fs]')
    plt.plot(oh, 'r-')
    
    plt.subplot(212)
    plt.xlabel('Number of bonds/O atom')
    plt.annotate("using t0 = {0} fs".format(t0_histo), xy=(0.75, 0.75),
                 xycoords="axes fraction")
    plt.hist(oh[t0_histo:], 10)
    if write:
        plt.savefig(md_type+'_oh_number.eps', format='eps')
        plt.savefig(pp, format='pdf')
    plt.figure(5)
    plt.title(title)
    plt.xlabel('O-H distances')
    plt.hist(oh_distances[t0_histo:],bin_histo)
    if write:
        plt.savefig(md_type+'_oh_distance.eps', format='eps')
        plt.savefig(pp, format='pdf')
    return 


verlet= True
berendsen = False
write = True

if verlet:
    temp = '500K'
    path_traj = './NVE/'+temp+'/'
    title = 'NVE Verlet 1nm PVA'
    md_type = 'verlet'
    t0_histo = 2000
    multipage = md_type+'_'+temp+'.pdf'

if berendsen:
    path_traj = './NVT/berendsen/'
    title = 'NVT berendsen 1nm PVA'
    md_type = 'berendsen'
    t0_histo = 2000
    multipage = md_type+'.pdf'
bin_histo = 20

pp = PdfPages(multipage)
plot_energies(path_traj,md_type,title,t0_histo, bin_histo,
              write=write)

bin_histo = 10
plot_oh(path_traj,md_type,title,t0_histo, bin_histo,
        write=write)

plt.show()
pp.close()
