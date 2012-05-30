#red blue green rbg, dashes squares triangles -- s ^
from ase.io import read
from ase.io.trajectory import PickleTrajectory
from ase.calculators.neighborlist import NeighborList
import matplotlib.pyplot as plt
from ase.units import fs, kB
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

def get_junction_vectors(positions,junction_index):
    l = junction_index
    l_f = l[1:]
    l_i = l[:-1]
    a_l = positions[:,l_f,:] - positions[:,l_i,:]
    return a_l

def plot_md(path,md_type,title,t0_histo,bin_histo,nf=1,nc=1,write=False):
    """ nc=1 by default caculates a_l for every carbon, 
            an integer every l carbons"""
    traj = PickleTrajectory(path+md_type+'.traj')
    #Get the C-C backbone position array.
    c_index = [ a.index for a in traj[0] if a.symbol == 'C']
    c_traj = [ atoms[c_index] for atoms in traj ]
    c_pos = np.array([ atoms.get_array('positions') for atoms in c_traj ])
    #define the junction points and get junction vectors
    print 'Number of C in system', len(c_traj[0])
    l = range(0,len(c_traj[0]),nc)
    a_l = get_junction_vectors(c_pos,l)
    #calculate end-to-end vector r_ee sum over junctions
    r_ee = a_l.sum(axis=1)
#    temp = [ atoms.get_temperature() for atoms in traj]
    r2 = (r_ee * r_ee).sum(axis=1)
    r2_mean = r2.mean()
    r = np.sqrt(r2)

    
    plt.figure(nf)
    plt.subplot(211)
    plt.title(title)
    plt.ylabel('R [Ang]')
    plt.xlabel('Time [fs]')
    plt.plot( r, 'g-' , label=r'R[t]')
    plt.legend(loc='upper right')

    plt.subplot(212)
    plt.xlabel('R [Ang]')
    plt.annotate("using t0 = {0} fs".format(t0_histo), xy=(0.75, 0.75),
                 xycoords="axes fraction")
    plt.hist(r[t0_histo:] , bin_histo)

    if write:
        #plt.savefig(md_type+'_R.eps', format='eps')
        plt.savefig(pp, format='pdf')

    return


verlet= True
berendsen =False
write = True

if verlet:
    temp = '500K' #300K init -> 260K mean, 500K init-> 320K mean
    path_traj = './NVE/'+temp+'/'
    title = 'NVE Verlet 1nm PVA 320K'
    md_type = 'verlet'
    t0_histo = 2000
    multipage = md_type+'_'+temp+'_R.pdf'
if berendsen:
    path_traj = './NVT/berendsen/'
    title = 'NVT berendsen 1nm PVA'
    md_type = 'berendsen'
    t0_histo = 2000
    multipage = md_type+'_R.pdf'

bin_histo = 20
nc_list = [2,4,6,8,10] #index increment = monomers * carbon_per_monomer 
nf = 1
pp = PdfPages(multipage)
for nc in nc_list:
    title_nc = title+' nc '+str(nc)
    plot_md(path_traj,md_type,title_nc,t0_histo, bin_histo,
            nf=nf,nc=nc,write=write)
    nf += 1


plt.show()

pp.close()
