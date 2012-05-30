from numpy import savetxt,shape,zeros,asanyarray,transpose,add,arange
from numpy import loadtxt,indices,reshape,ones,sum,where,sin,cos,exp,pi
from ase import *
from ase.parallel import paropen
from ase.units import Bohr, Hartree
from gpaw import *
from gpaw.utilities.dos import raw_orbital_LDOS, fold
from gpaw.analyse.expandyl import ExpandYl 
#from gpaw.utilities.tools import coordinates


class MPC:

    """ Class to perform current standard analysis of monolayer protecte clusters
    #init:
    calc = Calculator('out.gpw')
    mpc = MPC(calc)
    mpc.set_core('Au',coordination=[(n_symbol,n_max,d_max)])     #define core
    #use in serial:
    mpc.write_local_dos(list) #write sp-df dos for core(default) or list=[i,...]
    with i,... atom index
    mpc.get_radial_distribution(dR=dR)
    mpc.calculate_sp_band(plt=True,ylm=True) #set plt to write out plt
    #wavefunctions and ylm to expand the sp-band to find jellium states
    #default False

    #use in parallel (dos function work in serial only):
    mpc.read_summary()
    mpc.calculate_ylm_from_list(list=mpc.sp_band,Rmax=2.0)
    """

    def __init__(self, calc):
        
        self.gd = calc.density.gd
        self.atoms= calc.get_atoms()
        self.spinpol = calc.get_spin_polarized()
        self.nspins = calc.get_number_of_spins()
        self.density = calc.density
        self.energies = []
        self.energies += [calc.get_eigenvalues(spin=0)]
        if self.spinpol:
            self.energies += [calc.get_eigenvalues(spin=1)]
        self.calc = calc
        self.cm = self.atoms.get_center_of_mass()
        self.sp_band = None
        c = calc.input_parameters['convergence']['bands']
        if c == 'all' or c == 'occupied':
            c = len(self.energies[0])
        if c < 0 :
            c += len(self.energies[0])
        self.converged = c
        self.a0 = Bohr
        self.f_n = []
        self.f_n += [self.calc.get_occupation_numbers(kpt=0,spin=0)]
        if self.spinpol:
            self.f_n += [self.calc.get_occupation_numbers(kpt=0,spin=1)]
        self.w = None

    def set_core(self,symbol,coordination=None):
        self.core = atom_symbol(self.atoms,symbol,coordination)
        
    def write_local_dos(self,list=None,npts=5000,width=0.05,header='core'):
        angular=['s','p','d','f']
        if list is None:
            list = self.core
        e_shape = shape(self.energies)
        self.w = zeros((e_shape[0],len(angular),e_shape[1]))
        
        for s in range(self.nspins):
            out_name = 'dos_'+header+'_s'+str(s)+'.dat'
            out = paropen(out_name,mode='w')
            out.write('#'+header+' atoms, npts %d, width %f\n'%(npts,width))
            out.write('#energies eV \t angular=s \t p \t d \t f\n')
            for j,l in enumerate(angular):
                for i in list:
                    temp, temp2 = raw_orbital_LDOS(self.calc, i, s,l)
                    self.w[s,j,...] += temp2

            w_gauss = []
            for j,l in enumerate(angular):
                temp3,temp4 = fold(self.energies[s], 
                                   self.w[s,j,...], npts, width)
                if j == 0:
                    w_gauss = [temp3]
                w_gauss += [temp4]
            data = asanyarray(w_gauss)
            savetxt(out,transpose(data),delimiter='\t')

            
    def calculate_sp_band(self,wmin=50.0,w0=None,ylm=False,Rmax=None,plt=False):
        """ returns list w>wmin as % of wtot,
        write ylm analysis and write plt if set to true"""
        if w0 is None:
            w0 = 0.2 / self.nspins

        if self.sp_band is None:
            self.sp_band = []
            for s in range(self.nspins):
                list = []
                sp = add(self.w[s,0,...],self.w[s,1,...])
                df = add(self.w[s,2,...],self.w[s,3,...])
                tot = sp+df
                percentage = sp/tot*100.0
                for i in range(self.converged):
                    if (tot[i]>w0 and percentage[i]>wmin):
                        list += [i]
                self.sp_band += [list] 

        if plt == True:
            write('coords-plt.xyz',self.atoms)
            for s in range(self.nspins):
                write_wf(self.atoms,self.calc,self.sp_band[s],s)
        if ylm == True:
            self.calculate_ylm_from_list(self.sp_band,Rmax=Rmax)

        #write summary
        for s in range(self.nspins):
            g = paropen('summary'+str(s)+'.dat',mode='w')
            g.write('#spin\tband\tenergy[eV]\tocc\tsp_band(yes=1)\n')
            g.write('#wmin,w0= %f, %f\n' %(wmin,w0))
            g.write('#number of  orbitals in sp band:%d\n'%len(self.sp_band[s]))

            e = shape(self.energies)
            sp = [0.0]*e[1]
            for i in self.sp_band[s]:
                sp[i] =1.0
            data = [ [[0.0]*e[1],arange(e[1]),self.energies[0],
                      self.f_n[0],sp] ]
            data_out = asanyarray(data)
            help(savetxt)
            savetxt(g,
                    transpose(data_out),
                    delimiter='\t',
                    fmt='%8.4f ' )


    def read_summary(self):
        self.sp_band = []
        for s in range(self.nspins):
            g = paropen('summary'+str(s)+'.dat',mode='r')
            data_out = loadtxt(g,delimiter='\t',comments='#')
            data = transpose(data_out)
            list = []
            e = shape(data)
            for i in range(e[1]):
                if data[4,i] == 1.0:
                    list += [i]
            self.sp_band += [list]

    def write_bader(self):
        n = self.calc.get_all_electron_density(gridrefinement=2)
        n *= Bohr**3
        
        write('density.cube', self.atoms, data=n)
        

    def calculate_ylm_from_list(self,list=None,Rmax=None,atom_origin='cm'):
        self.calc.wfs.initialize_wave_functions_from_restart_file()
        f = paropen('ylm.dat',mode='w')
        if atom_origin == 'cm':
            ylm = ExpandYl(self.cm,self.gd,lmax=6,Rmax=Rmax)
        else:
            center = self.atoms.get_positions()[atom_origin]
            ylm = ExpandYl(center,self.gd,lmax=6,Rmax=Rmax)
        
        f.write("#YL expansion")
        f.write("#center of mass, Rmax=%f Ang, kpt fixed to 0\n"%Rmax)
        f.write("#spin\tband\t E(eV)\tocc\tnorm\tsum\tweight\t percentage(respect to weight)s-p-d-f-g-h-i \n")
        krange = range(len(self.calc.wfs.ibzk_kc))
        for s in range(self.nspins):
            for n in list[s]:
                for k in krange:
                    u = k*self.calc.wfs.nspins + s
                    kpt = self.calc.wfs.kpt_u[u]               
                    temp =  kpt.psit_nG[n]
                    norm = self.gd.integrate(temp**2)
                    gl,weight = ylm.expand(temp)

                    gsum = sum(gl)
                    gl = 100 * gl / weight
                    out = '%2d  %5d' % (s, n)
                    out += '%10.4f %8.4f' % (self.energies[s][n],
                                             self.f_n[s][n])
                    out += "%8.4f %8.4f %8.4f" % (norm, gsum, weight)
                    
                    for g in gl:
                        out += "%8.2f" %g
                    out += "\n"
                    f.write(out)
    def calculate_ylm_from_list_old(self,list=None,Rmax=None):
        ylm = ExpandYl(self.cm,self.calc.gd,lmax=6,Rmax=Rmax)
        ylm.to_file(self.calc,bands=list[0])

    def get_vacuum(self,s=0,delta=0.5):
        """ returns value of effective ks potential on the boundaries """
        #get effective potential
        v = self.calc.get_effective_potential(spin=s, pad=False)
        #get coordinates
        xyz, r2 = coordinates(self.gd,self.a0)
        #get slices:planes xy,xz,zy
        cell = self.atoms.get_cell()
        for i in range(3):
            plane = v[where(xyz[i,:,:,:]<=delta)]
        #get averages and variances
            nlen = len(plane)
        #print
            print nlen
            print 'mean:',plane.mean(),'st. dev.:',sqrt(plane.var())
        for i in range(3):
            plane = v[where(xyz[i,:,:,:]>=(cell[i,i]-delta))]
        #get averages and variances
            nlen = len(plane)
        #print
            print nlen
            print 'mean:',plane.mean(),'st. dev.:',sqrt(plane.var())
        return

    def get_radial_distribution(self,dR=1.0,save=True,mode='ae'):
        """ write rho(R) the radial all electron distribution """
        """ return bins left-right, histo[rho(r)]dr """
        dh = self.calc.get_grid_spacings()
        dv = dh[0]*dh[1]*dh[2]
        #get radial grid
        v = self.cm
        #get all electron density
        if mode == 'ae':
            rho = self.calc.get_all_electron_density(gridrefinement=2,
                                                     pad=False)
            dv /= 2**3
            xyz, r2 = coordinates(self.calc.finegd,self.calc.a0,v=v)

        elif mode == 'valence':
            rho = self.calc.get_pseudo_density(pad=False)
            xyz, r2 = coordinates(self.gd,self.calc.a0,v=v)
        #check:
        #print sum(rho)*dv
        
        rmax = round(sqrt(amax(r2)))+1.
        bin = (arange(0,rmax,dR))**2
        new = rho*dv
        newbin = bin[newaxis,:]
        histo,rarray = histogramdd(r2.ravel(),newbin,weights=new.ravel())
        if save == True:
            out = 'radial.dat'
            save_radial_file(out,histo,rarray[0],dR)
        return histo,rarray

    
def save_radial_file(out,histo,rarray,dR): #save to file  
    out = paropen(out,mode='w')
    out.write('# all electron density radial distribution\n')
    out.write('# first column:left-right bins spaced by dR %f Ang, rho(r)dr\n'%dR)
    out.write('# second: 4pi r**2 rho(r)dr -> sum of third column gives N total number of electrons\n')
    out.write('# third column:center of bins = r\n')
    out.write('# fourth c: rho(r)\n')
    out.write('# in gnuplot use: plot "radial.dat" u 1:4 w steps \n')
    data = zeros((4,len(rarray)))
    data[0,...] = np.sqrt(rarray)
    data[1,:-1] = histo
    data[2,...] = data[0,...]+dR/2.
    data[3,:-1] = histo/(4.*pi*data[2,:-1]**2)
    savetxt(out,transpose(data),delimiter='\t')


def atom_symbol(atoms,symbol,coordination=None):
    """ Example atom_symbol(cl,'Au',coordination=[('S',1,2.6)]) 
    returns core for most of MPC, Au102 needs one extra Au in center"""
    list = []
    for i,a in enumerate(atoms):
        if a.symbol == symbol:
            if coordination is None:
                list += [i]
            else:
                for n_symbol,n_max,d_max in coordination:
                    c = 0
                    for j,b in enumerate(atoms):
                        if b.symbol == n_symbol:
                            if atoms.get_distance(i,j)<=d_max:
                                c += 1
                    if c<=n_max:
                        list += [i]
                
    return list

def write_wf(atoms,calc,list,spin):
    for n in list:
        temp = calc.get_pseudo_wave_function(band=n,spin=spin)
        name = 'wf_%d_s%d.plt' % (n,spin)
        write(name,atoms,data=temp)


def coordinates(gd,a0,v=None):
    """Constructs and returns matrices containing cartesian coordinates,
       and the square of the distance from vector v
       If v=None the corner (0,0,0) of the box is taken.
       The origin is placed in the center of the box described by the given
       grid-descriptor 'gd'.
    """    
    I  = indices(gd.n_c)
    dr = reshape(gd.h_c, (3, 1, 1, 1))
    if v is None:
        #r0 = reshape(gd.h_c * gd.beg_c - .5 * gd.cell_c, (3, 1, 1, 1))
        r0 = reshape(gd.h_c * gd.beg_c , (3, 1, 1, 1))
    else:
        r0 = reshape(gd.h_c * gd.beg_c - v/a0, (3, 1, 1, 1))
    r0 = ones(I.shape) * r0
    xyz = r0 + I * dr
    xyz = xyz*a0
    r2 = sum(xyz**2, axis=0)


    # Remove singularity at origin and replace with small number
    #middle = gd.N_c / 2.
    # Check that middle is a gridpoint and that it is on this CPU
    #if (npy.alltrue(middle == npy.floor(middle)) and
    #    npy.alltrue(gd.beg_c <= middle) and
    #    npy.alltrue(middle < gd.end_c)):
    #    m = (middle - gd.beg_c).astype(int)
    #    r2[m[0], m[1], m[2]] = 1e-12

    # Return r^2 matrix
    return xyz, r2


def sum_optical(list):
    axis = list[0]
    g = paropen('spectrum_'+axis+'.dat',mode='r')
    data_out = loadtxt(g,comments='#')
    g.close()
    data = transpose(data_out)

    for i in range(1,len(list)):
        axis = list[i]
        g = paropen('spectrum_'+axis+'.dat',mode='r')
        temp_out = loadtxt(g,comments='#')
        g.close()
        temp = transpose(temp_out)

        for j in [1,2,3]:
            data[j,...] += temp[j,...]
            

    out = paropen('spectrum_sum.dat',mode='w')
    savetxt(out,transpose(data),delimiter='\t')

    
def write_density_difference(list1,list2):
    name1 = list1[0]
    name2 = list2[0]
    calc1 = Calculator(name1)
    atoms1 = calc1.get_atoms()
    calc2 = Calculator(name2)
    atoms2 = calc2.get_atoms()
    rho1 = calc1.get_all_electron_density(gridrefinement=2)
    rho2 = calc2.get_all_electron_density(gridrefinement=2)
    for name in list1[1:]:
        calc = Calculator(name)
        temp = calc.get_all_electron_density(gridrefinement=2)
        rho1 += temp
    for name in list2[1:]:
        calc = Calculator(name)
        temp = calc.get_all_electron_density(gridrefinement=2)
        rho2 += temp
        
    diff = rho2-rho1
    
    write('rho1.plt',atoms1,data=rho1)
    write('rho2.plt',atoms2,data=rho2)
    write('diff.plt',atoms2,data=diff)

def xrd(t,atoms,dlist,f,B,l,geometry):
    s = 2. * sin(t) / l
    n = len(atoms)
    sum = 0.0
    for i in range(n):
        for j in range(i+1,n):
            d = dlist[i,j]
            sin2pids = sin(2.*pi*d*s)
            fij = f[i]*f[j]
            sum += 2.* fij * sin2pids / (2.*pi*d*s)
    if B:
        expterm = exp(-B/2.*s**2)
        sum *= expterm
    if geometry:
        alpha = 1.01
        cost = cos(t)
        cos2t = cos(2*t)
        geo = cost / (1.+alpha*cos2t**2)
        sum *= geo
    return s,sum

def write_xrd(atoms_all,element='all',scale=True,B=False,l=0.15405,geometry=False,name='xrd.dat',bader=False):
    steps = 500
    o = paropen(name,'w')
    if element=='all':
        atoms = range(len(atoms_all))
    elif element=='core':
        atoms = atom_symbol(atoms_all,'Au',coordination=[('S',1,2.7)])
    else:
        atoms = atom_symbol(atoms_all,element)
    n = len(atoms)
    #print "atoms",len(atoms)
    f = atoms_all.get_atomic_numbers()
    if bader:
        bader_file = 'ACF.dat'
        attach_charges(atoms_all, bader_file,shift=True)
        c = atoms_all.get_charges()
        f = f*1.0 -c
    dlist = zeros([n,n])
    if scale:
        factor = 0.974
        #factor of 10 for distance in nm, factor 0.974 to account theory->exp bulk au distances
    else: 
        factor = 1.0

    for i in range(n):
        for j in range(i+1,n):
            dlist[i,j] =  atoms_all.get_distance(i,j)/10.*factor


    pre_factor = (f*f).sum()
    o.write('#t,s=2sin(t)/lambda,I(s)\n')
    for i in range(1,steps+1):
        t = i * pi/(2.*steps)
        s,sum = xrd(t,atoms,dlist,f,B,l,geometry)
        out = "%f\t%f\t%f\n" % (t,s,pre_factor+sum)
        o.write(out)

def delete(atoms,del_list):
    del_list.sort(reverse=True)
    for i in del_list:
        atoms.pop(i)

def write_list(atoms,list,file):
    new = Atoms()
    for i in list:
        new.append(atoms[i])
    write(file,new)

def create_pairs(lr,excitations,v_min,fe_sp):
    """ return the pairs having v>v_min in the given range"""
    new_list = []
    for i in excitations:
        ex = lr[i]
        #print ex.analyse()
        for f,k in zip(ex.f,ex.kss):
            if f*f > v_min:
                if k.i <= fe_sp:
                    #mark sp with -1
                    a = -1
                else:
                    a = k.i
                if k.j <= fe_sp:
                    #mark sp with -1
                    b = -1
                else:
                    b = k.j
                pair = (a,b)
                new_list.append(pair)
    reduced = set(new_list)
    return sorted(list(reduced))


def get_osc(lr):
    """Get excitation energies in Hartrees"""
    el = []
    for ex in lr:
        el.append(ex.GetOscillatorStrength()[0])
    return array(el)
    
def GetOscillator(lr):
    """Get excitation energies in units of the calculators energy unit"""
    return get_osc(lr) * 1.0

def range_excitations(emax,lr,osc_min):
    energies = lr.GetEnergies()
    excs = where(energies<emax)[0]
    osc = GetOscillator(lr)
    list = []
    for ex in excs:
        if osc[ex]>osc_min:
            list.append(ex)
    print "numer of excitations",len(list)
    return list

def write_pairs(file,v_min,e_max,fe_sp,pairs,excitations,lr):
    """ write"""
    g = open(file,'w')
    g.write("#Write the principal eigenvector weights") 
    g.write("(rest sum to 1) with f*f> %.5f \n" %v_min)
    g.write("#e_max:%.3f eV\n"%e_max)
    g.write("#energy eV, pairs (-1=Fe sp band wf<%d)\n"%fe_sp)
    i=2
    g.write("#")
    for p in pairs:
        g.write(str(i)+":"+str(p)+" ")
        i +=1
    g.write("\n")
    lp = len(pairs)            
    out = zeros(lp)
    for i in excitations:
        ex = lr[i]
        g.write('%.3f'%(ex.energy*27.211))
        for f,k in zip(ex.f,ex.kss):
            if f*f > v_min:
                if k.i <= fe_sp:
                    #mark sp with -1
                    a = -1
                else:
                    a = k.i
                if k.j <= fe_sp:
                    #mark sp with -1
                    b = -1
                else:
                    b = k.j
                pair = (a,b)
                p = pairs.index(pair)
                for m in range(p,lp):
                    out[m] += f*f
        for a in out:
            g.write("\t%.3f"%a)
        g.write("\n")
        out[:] = 0.0
    ll = range(2,lp+2)
    for i in ll[::-1]:
        print "'v_sum.dat' u 1:%d  w i  t '%d'  ls %d,\ " %(i,i,i)
    
def write_excitations(lr=None,e_max=1.6,file='v_sum.dat',
                      osc_min=0.005,v_min=0.05,fe_sp=-1):
    excitations = range_excitations(e_max,lr,osc_min) 
    pairs = create_pairs(lr,excitations,v_min,fe_sp)
    write_pairs(file,v_min,e_max,fe_sp,pairs,excitations,lr)
    return excitations

def reduce(l):
    reduced = set(l)
    return sorted(list(reduced))

def get_pairs_list(list1,list2):
    """ returns the pairs (i,j) excluding cases (i,i) """
    l = []
    for i in list1:
        for j in list2:
            if i != j:
                pair = sorted([i,j])
                l.append((pair[0],pair[1]))
    return reduce(l)

def get_distance_array(list,atoms):
    n = len(list)
    array = zeros(n)
    for i in range(n):
        x = list[i]
        array[i] = atoms.get_distance(x[0],x[1])
    return array

def get_mean_var(array,bmax='all',rmax=None):
    if bmax=='all':
        bmax=len(array)

    if rmax != None:
        narray = array[where(array<rmax)]
        b = len(narray)
    else:
        b = 0
        c = 1.0
        while b < bmax:
            c += 0.0001
            narray = array[where(array<c)]
            b = len(narray)
    print b,bmax,narray.mean(),sqrt(narray.var())

def  write_array(array,width,filename,npts=200,weights=None,):
    from gpaw.utilities.dos import fold
    if weights==None:
        weights  = np.ones(len(array))
    x,y = fold(array, weights, npts, width)
    out = open(filename,'w')
    savetxt(out,zip(x,y),delimiter='\t')

def sum_conductance(list):
    k = str(list[0])
    g = paropen('k'+k+'/T_0'+k+'.dat',mode='r')
    data_out = loadtxt(g,comments='#')
    g.close()
    data = transpose(data_out)

    for i in range(1,len(list)):
        k = str(list[i])
        g = paropen('k'+k+'/T_0'+k+'.dat',mode='r')
        temp_out = loadtxt(g,comments='#')
        g.close()
        temp = transpose(temp_out)        
        data[1,:] += temp[1,:]
            
    data[1,:] /= len(list)

    out = paropen('conductance.dat',mode='w')
    savetxt(out,transpose(data),delimiter='\t')


def print_bader(atom_file='coords-plt.xyz',fileobj='ACF.dat'):
    from ase.io.bader import attach_charges
    atoms = read(atom_file)
    attach_charges(atoms, fileobj)
    g = paropen('bader.dat','w')
    out = 'Atom \t Bader charge \t Radius Ang\n'
    g.write(out)
    cm = atoms.get_center_of_mass()
    for atom in atoms:
        r2 = ((atom.get_position()-cm)**2).sum()
        out = '%s \t %.2f \t %.2f \n' % (atom.symbol, atom.charge,sqrt(r2))
        print out
        g.write(out)

def attach_charges_mine(atoms, fileobj='ACF.dat', displacement=1e-4,shift=False):
    import numpy as np
    from ase.units import Bohr
    """Attach the charges from the fileobj to the Atoms."""
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    sep = '---------------'
    while sep not in fileobj.readline(): # Skip to after first seperator line                                               
        pass

    for line in fileobj:
        if sep in line: # Stop at last seperator line                                                                       
            break

        words = line.split()
        if len(words) != 6:
            raise IOError('Number of columns in ACF file incorrect!\n'
                          'Check that Bader program version >= 0.25')

        atom = atoms[int(words[0]) - 1]
        atom.charge = atom.number - float(words[4])

        if displacement is not None: # check if the atom positions match                                                    
            xyz = np.array([float(w) for w in words[1:4]]) * Bohr
            assert np.linalg.norm(atom.position - xyz) < displacement
        if shift:
            xyz = np.array([float(w) for w in words[1:4]]) * Bohr
            atom.set_position(xyz)
