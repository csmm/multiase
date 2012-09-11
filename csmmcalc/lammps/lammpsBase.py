#!/usr/bin/env python

# lammps.py (2011/03/29)
# An ASE calculator for the LAMMPS classical MD code available from
#       http://lammps.sandia.gov/
# The environment variable LAMMPS_COMMAND must be defined to point to the LAMMPS binary.
#
# Copyright (C) 2009 - 2011 Joerg Meyer, joerg.meyer@ch.tum.de
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this file; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# or see <http://www.gnu.org/licenses/>.


import os 
import shutil
import shlex
import time
from subprocess import Popen, PIPE
from threading import Thread
from tempfile import mkdtemp, NamedTemporaryFile
import numpy as np
import decimal as dec
from ase import Atoms
import lammpsUnitConversion

__ALL__ = ['LAMMPSBase']

# "End mark" used to indicate that the calculation is done
CALCULATION_END_MARK = '__end_of_ase_invoked_calculation__'

class LAMMPSParameters:
    def __init__(self, **pars):
        self.atom_style     = pars.get('atom_style', 'full')
        self.units          = pars.get('units', 'metal')
        self.neighbor       = pars.get('neighbor')
        self.newton         = pars.get('newton')
        
        self.pair_style     = pars.get('pair_style')
        self.bond_style     = pars.get('bond_style')
        self.angle_style    = pars.get('angle_style')        
        self.dihedral_style = pars.get('dihedral_style')
        self.improper_style = pars.get('improper_style')
        self.special_bonds  = pars.get('special_bonds')
        
        # Pair coeffs to be specified in the input file
        self.pair_coeffs    = pars.get('pair_coeffs', [])
        
        self.extra_cmds     = pars.get('extra_cmds', [])

class LAMMPSData:
    def __init__(self):
        self.clear()
        
    def clear(self):
        self.masses = None
        self.coeff_tables = []
        self.atom_types = None
        self.bonds = None
        self.angles = None
        self.dihedrals = None
        self.impropers = None
        
    def add_coeffs(self, title, table):
        self.coeff_tables.insert(dict(title=title, table=table))
        
        
        
class LAMMPSBase:

    def __init__(self, label='lammps', tmp_dir=None, parameters={}, 
                 update_charges = False,
                 keep_alive=True, debug=False):
        """The LAMMPS calculators object """
  
        self.label = label
        self.parameters = LAMMPSParameters(**parameters)
        self.data = LAMMPSData()
        self.calls = 0
        self.forces = None
        self.atoms = None
        self.error = None
        self.update_charges = update_charges
        self.keep_alive = keep_alive
        self.keep_tmp_files = debug           
        self._lmp_handle = None        # To handle the lmp process

        self._custom_thermo_args = ['step', 'temp', 'press', 'cpu', 
                                    'pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz',
                                    'ke', 'pe', 'etotal',
                                    'vol', 'lx', 'ly', 'lz', 'atoms']
        
        # Thermo output of LAMMPS, a list of dictionaries. See read_lammps_output()
        self.thermo_content = []
        
        self._dump_fields = ['id', 'type', 'x', 'y', 'z', 'vx',
                             'vy', 'vz', 'fx', 'fy', 'fz', 'q']
        
        if tmp_dir is None:
            self.tmp_dir = mkdtemp(prefix='LAMMPS-')
        else:
            # If tmp_dir is pointing somewhere, don't remove stuff!
            self.keep_tmp_files = True
            self.tmp_dir=os.path.realpath(tmp_dir)
            if not os.path.isdir(self.tmp_dir):
                os.mkdir(self.tmp_dir, 0755)
        
        if debug:
           print 'LAMMPS (label: %s) running at %s' % (self.label, self.tmp_dir)
     
    def prepare_calculation(self):
        """ Implement this method in subclasses """
        raise NotImplementedException()
        
        
    def use_ff_file(self, filepath, target_filename=None):
        if not target_filename: 
            target_filename=os.path.split(filepath)[1]
        shutil.copy(filepath, os.path.join(self.tmp_dir, target_filename))

    def clean(self, force=False):
        self._lmp_end()
        self.data.clear()
        # TODO: When to delete tmp dir?
        #if not self.keep_tmp_files:
        #    shutil.rmtree(self.tmp_dir)

    def get_potential_energy(self, atoms):
        self.update(atoms)
        energy = self.thermo_content[-1]['pe']
        return self.to_ase_units(energy, 'energy')

    def get_forces(self, atoms):
        self.update(atoms)
        return self.to_ase_units(self.forces, 'force')

    def get_stress(self, atoms):
        self.update(atoms)
        tc = self.thermo_content[-1]
        stress =  np.array([tc[i] for i in ('pxx','pyy','pzz',
                                         'pyz','pxz','pxy')])
        return self.to_ase_units(stress, 'stress')  
    
    def calculation_required(self, atoms, quantities=None):
        return atoms != self.atoms
        
    def update(self, atoms):
        if not self.calculation_required(atoms): return
        self.setup_calculation(atoms)
        self.run_single_step()
        self.read_lammps_output()
        self.read_lammps_trj()
        self.close_calculation()
    
    def minimize(self, atoms, etol=0, ftol=0, maxeval=100000, maxiter=100000000, min_style=None):
        if etol == 0 and ftol == 0: 
            raise RuntimeError('Specify at least one tolerance value!')
        self.update(atoms)
        ftol = self.from_ase_units(ftol, 'force')
        minimize_params = '%g %g %i %i' % (etol, ftol, maxeval, maxiter)
        self.atoms = atoms
        self.setup_calculation(atoms)
        self.run_minimization(minimize_params, min_style)
        self.read_lammps_output()
        self.read_lammps_trj()
        self.close_calculation()

        
    def molecular_dynamics(self, atoms, timestep, fix):
        self.update(atoms)
        fix = fix
        timestep = str(self.from_ase_units(timestep, 'time'))        
        self.setup_calculation(atoms)
        nsteps = yield
        cur_step = 0
        while nsteps:
            self.set_dumpfreq(cur_step+nsteps)
            cur_step += nsteps
            self.run_md(fix, timestep, nsteps)
            self.read_lammps_output()
            self.read_lammps_trj()
            nsteps = yield
        self.close_calculation()
        
    def setup_calculation(self, atoms):
        inv_cell = np.linalg.inv(atoms.cell)
        frac_positions = np.dot(inv_cell, atoms.positions.T)
        if np.any(frac_positions < 0) or np.any(frac_positions > 1):
	  raise RuntimeError("Atoms not inside calculation cell! Use Atoms.center(vacuum).")
	
        self.atoms = atoms
        self.prepare_calculation()
        self.prism = prism(self.atoms.get_cell())
        self.prepare_lammps_run()
    
    def close_calculation(self):
        self.finish_lammps_run()
        self.clean()
        self.atoms = self.atoms.copy()
        
    def prepare_lammps_run(self):
        """ Start LAMMPS and prepare everything for 'run' or 'minimize' commands """
        if not self._lmp_alive(): self.start_lammps_process()
        self.prepare_lammps_io()
        self.write_lammps_input()
        
        
    def start_lammps_process(self):
        # set LAMMPS command from environment variable
        if 'LAMMPS_COMMAND' in os.environ:
            lammps_cmd_line = shlex.split(os.environ['LAMMPS_COMMAND'])
            if len(lammps_cmd_line) == 0:
                self.clean()
                raise RuntimeError('The LAMMPS_COMMAND environment variable '
                                   'must not be empty')
            # want always an absolute path to LAMMPS binary when calling from self.dir                       
            lammps_cmd_line[0] = os.path.abspath(lammps_cmd_line[0])
        else:
            self.clean()
            raise RuntimeError('Please set LAMMPS_COMMAND environment variable')
        if 'LAMMPS_OPTIONS' in os.environ:
            lammps_options = shlex.split(os.environ['LAMMPS_OPTIONS'])
        else:
            lammps_options = shlex.split('-echo log -screen none')
            
        self._lmp_handle = Popen(lammps_cmd_line+lammps_options+['-log', '/dev/stdout'], 
                                    cwd=self.tmp_dir, stdin=PIPE, stdout=PIPE)
           
           
    def prepare_lammps_io(self):
        # setup input and output files for LAMMPS
        label = '%s%06d' % (self.label, self.calls)
        
        self.lammps_trj_file = NamedTemporaryFile(mode='r', prefix='trj_'+label,
                                    dir=self.tmp_dir, delete=(not self.keep_tmp_files))
        
        self.lammps_inputdata_file = NamedTemporaryFile(prefix='data_'+label,
                                    dir=self.tmp_dir, delete=(not self.keep_tmp_files))

        if self.keep_tmp_files:
            # Save LAMMPS input and output for reference
            inlog = NamedTemporaryFile(prefix='in_'+label, dir=self.tmp_dir, delete=False)
            outlog = NamedTemporaryFile(prefix='log_'+label, dir=self.tmp_dir, delete=False)
            
            self.lammps_input = special_tee(self._lmp_handle.stdin, inlog)
            self.lammps_output = special_tee(self._lmp_handle.stdout, outlog)
            
        else:
            self.lammps_output = self._lmp_handle.stdout
            self.lammps_input = self._lmp_handle.stdin
            
    
    def run_single_step(self):
        self.set_dumpfreq(1)
        self.lammps_input.write('run 0\n')
        self.lammps_input.write('print "%s"\n' % CALCULATION_END_MARK)
        self.lammps_input.write('log /dev/stdout\n') # Force LAMMPS to flush log
        
    def run_minimization(self, param_string, min_style=None):
        self.set_dumpfreq(999999)
        if min_style:
            self.lammps_input.write('min_style %s\n' % min_style)
        self.lammps_input.write('minimize %s\n'  % param_string)
        self.lammps_input.write('print "%s"\n' % CALCULATION_END_MARK)
        self.lammps_input.write('log /dev/stdout\n') # Force LAMMPS to flush log
        
    def run_md(self, fix, timestep, nsteps):
        self.lammps_input.write('fix fix_all %s\n' % fix)
        self.lammps_input.write('timestep %s\n' % timestep)
        self.lammps_input.write('run %s\n' % nsteps)
        self.lammps_input.write('print "%s"\n' % CALCULATION_END_MARK)
        self.lammps_input.write('log /dev/stdout\n') # Force LAMMPS to flush log
            
    def finish_lammps_run(self):
        if not self.keep_alive:
            self._lmp_end()

        exitcode = self._lmp_handle.poll()
        if exitcode and exitcode != 0:
            raise RuntimeError('LAMMPS exited in %s with exit code: %d.' %\
                                   (self.tmp_dir, exitcode))

        self.lammps_trj_file.close()
        self.lammps_inputdata_file.close()
        if self.keep_tmp_files:
            self.lammps_input.close_logfile()
            self.lammps_output.close_logfile()

            
    def write_lammps_input(self):
        """Write LAMMPS parameters and input data """
        f = self.lammps_input
        parameters = self.parameters
        
        f.write('# (written by ASE)\n')
        f.write('atom_style %s \n' % parameters.atom_style)
        
        pbc = self.atoms.get_pbc()
        f.write('units %s \n' % parameters.units)
        f.write('boundary %c %c %c \n' % tuple('sp'[x] for x in pbc))
        if parameters.neighbor:
            f.write('neighbor %s \n' % (parameters.neighbor))
        if parameters.newton:
            f.write('newton %s \n' % (parameters.newton))

        # Write interaction stuff
        f.write('\n### interactions \n')
        if parameters.pair_style:
            f.write('pair_style %s \n' % parameters.pair_style)
            
        if parameters.bond_style and len(self.data.bonds) != 0:
            f.write('bond_style %s \n' % parameters.bond_style)
            
        if parameters.angle_style and len(self.data.angles) != 0:
            f.write('angle_style %s \n' % parameters.angle_style)
            
        if parameters.dihedral_style and len(self.data.dihedrals) != 0:
            f.write('dihedral_style %s \n' % parameters.dihedral_style)
            
        if parameters.improper_style and len(self.data.impropers) != 0:
            f.write('improper_style %s \n' % parameters.improper_style)
        
        if parameters.special_bonds:
            f.write('special_bonds %s \n' % parameters.special_bonds)
        
        self.write_lammps_data()
        f.write('\n### read data \n')
        f.write('read_data %s\n' % self.lammps_inputdata_file.name)
        
        # Extra pair coeffs
        for line in parameters.pair_coeffs:
            f.write('pair_coeff %s \n' % line)
        
        for cmd in parameters.extra_cmds:
            f.write(cmd + '\n')
        
        f.write(('thermo_style custom %s\n' +
                'thermo_modify flush yes\n' +
                'thermo 0\n') % (' '.join(self._custom_thermo_args)))
                
        f.write('\ndump dump_all all custom ' +
                '1 %s %s\n' % (self.lammps_trj_file.name, ' '.join(self._dump_fields)) )

    def set_dumpfreq(self, freq):
        self.lammps_input.write('dump_modify dump_all every %s\n' % freq)
        
        
    def write_lammps_data(self):
        """Write system configuration and force field parameters to file to be read
           with read_data by LAMMPS."""
        
        f = self.lammps_inputdata_file
        data = self.data
        
        f.write(f.name + ' (written by ASE) \n\n')

        f.write('%d \t atoms \n' % len(data.atom_types))
        if data.bonds:     f.write('%d \t bonds \n' % len(data.bonds))
        if data.angles:    f.write('%d \t angles \n' % len(data.angles))
        if data.dihedrals: f.write('%d \t dihedrals \n' % len(data.dihedrals))
        if data.impropers: f.write('%d \t impropers \n' % len(data.impropers))
        
        def create_set_of_types(lst):
            if not lst: return None
            return set(element[0] for element in lst)
        
        atom_types     = set(self.data.atom_types)
        bond_types     = create_set_of_types(data.bonds)
        angle_types    = create_set_of_types(data.angles)
        dihedral_types = create_set_of_types(data.dihedrals)
        improper_types = create_set_of_types(data.impropers)

        f.write('%d  atom types\n' % len(atom_types))
        if data.bonds:     f.write('%d  bond types\n' % len(bond_types))
        if data.angles:    f.write('%d  angle types\n' % len(angle_types))
        if data.dihedrals: f.write('%d  dihedral types\n' % len(dihedral_types))
        if data.impropers: f.write('%d  improper types\n' % len(improper_types))

        p = prism(self.atoms.get_cell())
        xhi, yhi, zhi, xy, xz, yz = p.get_lammps_prism()

        f.write('0.0 %f  xlo xhi\n' % xhi)
        f.write('0.0 %f  ylo yhi\n' % yhi)
        f.write('0.0 %f  zlo zhi\n' % zhi)
        
        if p.is_skewed():
            f.write('%f %f %f  xy xz yz\n' % (xy, xz, yz))
        f.write('\n\n')
     
        def print_table(title, table):
            if len(table) == 0: return
            f.write('%s \n\n' % title)
            for index, row in enumerate(table, 1):
                f.write(('%d'+' %s'*len(row) +'\n') % ((index,) + tuple(row)))
            f.write('\n\n')
        
        if data.masses:
            print_table('Masses', [[mass] for mass in data.masses])
        
        for coeff_table in self.data.coeff_tables:
            print_table(coeff_table['title'], coeff_table['data'])
        
        # Atoms
        charges = self.atoms.get_charges()
        lammps_positions = p.pos_to_lammps(self.atoms.positions)
        positions = self.from_ase_units(lammps_positions, 'distance')
        if self.parameters.atom_style == 'full':
            table = [['1', tp, q] + list(pos)
                for tp, q, pos in zip(self.data.atom_types, charges, positions)]
        elif  self.parameters.atom_style == 'charge':
            table = [[tp, q] + list(pos)
                for tp, q, pos in zip(self.data.atom_types, charges, positions)]
        else:
            raise RuntimeError('Unsupported atom_style: %s' % self.parameters.atom_style)
        print_table('Atoms', table)
        
        if data.bonds:
            print_table('Bonds', data.bonds)
        if data.angles:
            print_table('Angles', self.data.angles)
        if data.dihedrals:
            print_table('Dihedrals', self.data.dihedrals)
        if data.impropers:
            print_table('Impropers', self.data.impropers)
        
        if self.atoms.has('momenta'):
            vel = self.atoms.get_velocities()
            lammps_velocities = self.from_ase_units(vel, 'velocity')
            print_table('Velocities', lammps_velocities)
        
        f.flush()
        
    def read_lammps_output(self):
        """ Read thermo output from LAMMPS stdout """
        f = self.lammps_output
        
        thermo_mark = ' '.join([x.capitalize() for x in
                                             self._custom_thermo_args[0:3]])
        
        thermo_content = []
        line = f.readline()
        while line and line.strip() != CALCULATION_END_MARK:
            if 'ERROR:' in line:
                raise RuntimeError('LAMMPS execution in %s failed. LAMMPS %s' 
                        % (self.tmp_dir,line))
            if line.startswith(thermo_mark):
                # get thermo output
                while True:
                    line = f.readline()
                    fields = line.split()
                    if len(fields) != len(self._custom_thermo_args): break
                    try:
                        # create a dictionary between each of the thermo_style args
                        # and it's corresponding value
                        thermo_content.append(dict(zip(self._custom_thermo_args, 
                                                       map(float, fields))))
                    except ValueError:
                        # Wasn't a thermo line after all
                        break
            else:
                line = f.readline()
        self.thermo_content = thermo_content

        
    def read_lammps_trj(self):
        """Read the LAMMPS dump file."""
        f = self.lammps_trj_file
        
        while True:
            line = f.readline()
            if not line:
                break

            if 'ITEM: TIMESTEP' in line:
                timestep = int(f.readline())
                n_atoms = 0
                lo = [] ; hi = [] ; tilt = []
                id = [] ; type = []
                positions  = np.empty((len(self.atoms), 3)) * np.nan
                forces     = np.empty((len(self.atoms), 3)) * np.nan
                velocities = np.empty((len(self.atoms), 3)) * np.nan
                charges    = np.empty((len(self.atoms), )) * np.nan

            if 'ITEM: NUMBER OF ATOMS' in line:
                line = f.readline()
                n_atoms = int(line.split()[0])
            
            if 'ITEM: BOX BOUNDS' in line:
                # save labels behind "ITEM: BOX BOUNDS" in triclinic case (>=lammps-7Jul09)
                tilt_items = line.split()[3:]
                for i in range(3):
                    line = f.readline()
                    fields = line.split()
                    lo.append(float(fields[0]))
                    hi.append(float(fields[1]))
                    if (len(fields) >= 3):
                        tilt.append(float(fields[2]))
            
            if 'ITEM: ATOMS' in line:
                # (reliably) identify values by labels behind "ITEM: ATOMS" - requires >=lammps-7Jul09
                column_names = line.split()[2:]
                for n in range(n_atoms):
                    line = f.readline()
                    fields = dict(zip(column_names, line.split()))
                    id = int(fields['id'])
                    
                    pos   = [ float(fields[col]) for col in ['x', 'y', 'z'] ]
                    vel   = [ float(fields[col]) for col in ['vx', 'vy', 'vz'] ]
                    force = [ float(fields[col]) for col in ['fx', 'fy', 'fz'] ]
                    q     = float(fields['q'])
                                        
                    rotate = self.prism.pos_to_ase
                    positions[id-1]  = rotate(self.to_ase_units(np.array(pos), 'distance'))
                    velocities[id-1] = rotate(self.to_ase_units(np.array(vel), 'velocity'))
                    forces[id-1]     = rotate(self.to_ase_units(np.array(force), 'force'))
                    charges[id-1]    = self.to_ase_units(np.array(q), 'charge')
                               
                self.forces = forces
                if timestep > 0:
                    self.atoms.set_positions(positions)
                    self.atoms.set_velocities(velocities)
                
                if self.update_charges:
                    self.atoms.set_charges(charges)
       
    
    def to_ase_units(self, value, quantity):
        return lammpsUnitConversion.convert(value, quantity, self.parameters.units, 'ASE')

    def from_ase_units(self, value, quantity):
        return lammpsUnitConversion.convert(value, quantity, 'ASE', self.parameters.units)

    def _lmp_alive(self):
        # Return True if this calculator is currently handling a running lammps process
        return self._lmp_handle and self._lmp_handle.poll() == None

    def _lmp_end(self):
        # Close lammps input and wait for lammps to end. Return process return value
        if self._lmp_alive():
            self._lmp_handle.stdin.close()
            return self._lmp_handle.wait()
        
        
class special_tee:
    """ Used for logging input and output to files """
    
    def __init__(self, orig_fd, out_fd):
        self._orig_fd = orig_fd
        self._out_fd = out_fd
        self.name = orig_fd.name
    def write(self, data):
        self._orig_fd.write(data)
        self._out_fd.write(data)
    def read(self, *args, **kwargs):
        data = self._orig_fd.read(*args, **kwargs)
        self._out_fd.write(data)
        return data
    def readline(self, *args, **kwargs):
        data = self._orig_fd.readline(*args, **kwargs)
        self._out_fd.write(data)
        return data
    def readlines(self, *args, **kwargs):
        data = self._orig_fd.readlines(*args, **kwargs)
        self._out_fd.write(''.join(data))
        return data
    def flush(self):
        self._orig_fd.flush()
        self._out_fd.flush()
        
    def close(self):
        self._orig_fd.close()
        self._out_fd.close()
        
    def close_logfile(self):
        self._out_fd.close()


class prism:
    def __init__(self, cell, pbc=(True,True,True), digits=10):
        """Create a lammps-style triclinic prism object from a cell
        """
        a, b, c = cell
        an, bn, cn = [np.linalg.norm(v) for v in cell]
        
        alpha = np.arccos(np.dot(b, c)/(bn*cn))
        beta  = np.arccos(np.dot(a, c)/(an*cn))
        gamma = np.arccos(np.dot(a, b)/(an*bn))
        
        xhi = an
        xyp = np.cos(gamma)*bn
        yhi = np.sin(gamma)*bn
        xzp = np.cos(beta)*cn
        yzp = (bn*cn*np.cos(alpha) - xyp*xzp)/yhi
        zhi = np.sqrt(cn**2 - xzp**2 - yzp**2)
    
        # Set precision
        self.car_prec = dec.Decimal('10.0') ** \
            int(np.floor(np.log10(max((xhi,yhi,zhi))))-digits)
        self.dir_prec = dec.Decimal('10.0') ** (-digits)
        self.acc = float(self.car_prec)
        self.eps = np.finfo(xhi).eps

        # For rotating positions from ase to lammps
        Apre = np.array(((xhi, 0,   0),
                         (xyp, yhi, 0),
                         (xzp, yzp, zhi)))
        self.R = np.dot(np.linalg.inv(cell), Apre)
        self.Rinv = np.linalg.inv(self.R)

        # Actual lammps cell may be different from what is used to create R
        def fold(vec, pvec, i):
            p = pvec[i]
            x = vec[i] + 0.5*p
            n = (np.mod(x, p) - x)/p
            return [float(self.f2qdec(a)) for a in (vec + n*pvec)]

        Apre[1,:] = fold(Apre[1,:], Apre[0,:], 0)
        Apre[2,:] = fold(Apre[2,:], Apre[1,:], 1)
        Apre[2,:] = fold(Apre[2,:], Apre[0,:], 0)

        self.A = Apre
        self.Ainv = np.linalg.inv(self.A)

        if self.is_skewed() and \
                (not (pbc[0] and pbc[1] and pbc[2])):
            raise RuntimeError('Skewed lammps cells MUST have '
                               'PBC == True in all directions!')
  
    def f2qdec(self, f):
        return dec.Decimal(repr(f)).quantize(self.car_prec, dec.ROUND_DOWN)

    def get_lammps_prism(self):
        A = self.A
        return (A[0,0], A[1,1], A[2,2], A[1,0], A[2,0], A[2,1])

    def pos_to_lammps(self, position):
        return np.dot(position, self.R)
    
    def pos_to_ase(self, position):
        return np.dot(position, self.Rinv)
    
    def is_skewed(self):
        acc = self.acc
        prism = self.get_lammps_prism()
        axy, axz, ayz = [np.abs(x) for x in prism[3:]]
        return (axy >= acc) or (axz >= acc) or (ayz >= acc)
        