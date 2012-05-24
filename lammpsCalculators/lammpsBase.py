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
from re import compile as re_compile, IGNORECASE
from tempfile import mkdtemp, NamedTemporaryFile, mktemp as uns_mktemp
import numpy as np
import decimal as dec
from ase import Atoms
from ase.parallel import paropen
from ase.units import GPa
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
        
        # Pair coeffs to be specified in the input file
        self.pair_coeffs    = pars.get('pair_coeffs', [])
            
        self.minimize       = pars.get('minimize')
        self.run            = pars.get('run')

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
                 files=[], always_triclinic=False, 
                 keep_alive=True, keep_tmp_files=False):
        """The LAMMPS calculators object

        keep_tmp_files: bool
            Retain any temporary files created. Mostly useful for debugging.
        tmp_dir: str
            path/dirname (default None -> create automatically).
            Explicitly control where the calculator object should create 
            its files. Using this option implies 'keep_tmp_files'
        no_data_file: bool
            Controls whether an explicit data file will be used for feeding
            atom coordinates into lammps. Enable it to lessen the pressure on
            the (tmp) file system. THIS OPTION MIGHT BE UNRELIABLE FOR CERTIAN
            CORNER CASES (however, if it fails, you will notice...).
        keep_alive: bool
            When using LAMMPS as a spawned subprocess, keep the subprocess
            alive (but idling whn unused) along with the calculator object.
        always_triclinic: bool
            Force use of a triclinic cell in LAMMPS, even if the cell is
            a perfect parallelepiped.
        """

        self.label = label
        self.parameters = LAMMPSParameters(**parameters)
        self.data = LAMMPSData()
        self.files = files
        self.calls = 0
        self.forces = None
        self.keep_alive = keep_alive
        self.keep_tmp_files = keep_tmp_files
        if tmp_dir is not None:
            # If tmp_dir is pointing somewhere, don't remove stuff!
            self.keep_tmp_files = True
        self._lmp_handle = None        # To handle the lmp process

        # read_log depends on that the first (three) thermo_style custom args
        # can be capitilized and matched aginst the log output. I.e. 
        # don't use e.g. 'ke' or 'cpu' which are labeled KinEng and CPU.
        self._custom_thermo_args = ['step', 'temp', 'press', 'cpu', 
                                    'pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz',
                                    'ke', 'pe', 'etotal',
                                    'vol', 'lx', 'ly', 'lz', 'atoms']
        self._custom_thermo_mark = ' '.join([x.capitalize() for x in
                                             self._custom_thermo_args[0:3]])

        # Match something which can be converted to a float
        f_re = r'([+-]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?|nan|inf))'
        n = len(self._custom_thermo_args)
        # Create a re matching exactly N white space separated floatish things
        self._custom_thermo_re = re_compile(r'^\s*' + r'\s+'.join([f_re]*n) + r'\s*$',
                                            flags=IGNORECASE)
        # thermo_content contains data "writen by" thermo_style.
        # It is a list of dictionaries, each dict (one for each line
        # printed by thermo_style) contains a mapping between each 
        # custom_thermo_args-argument and the corresponding
        # value as printed by lammps. thermo_content will be 
        # re-populated by the read_log method.
        self.thermo_content = []
        
        if tmp_dir is None:
            self.tmp_dir = mkdtemp(prefix='LAMMPS-')
        else:
            self.tmp_dir=os.path.realpath(tmp_dir)
            if not os.path.isdir(self.tmp_dir):
                os.mkdir(self.tmp_dir, 0755)
        
        for f in files:
            shutil.copy(f, os.path.join(self.tmp_dir, f))

    def clean(self, force=False):
        self._lmp_end()
        self.data.clear()

        if not self.keep_tmp_files:
            shutil.rmtree(self.tmp_dir)

    def get_potential_energy(self, atoms):
        self.update(atoms)
        energy = self.thermo_content[-1]['pe']
        return lammpsUnitConversion.toASE(energy, 'energy', self.parameters.units)

    def get_forces(self, atoms):
        self.update(atoms)
        return lammpsUnitConversion.toASE(self.forces, 'force', self.parameters.units)

    def get_stress(self, atoms):
        self.update(atoms)
        tc = self.thermo_content[-1]
        stress =  np.array([tc[i] for i in ('pxx','pyy','pzz',
                                         'pyz','pxz','pxy')])
        lammpsUnitConversion.toASE(stress, 'stress', self.parameters.units)

    def update(self, atoms):
        if not hasattr(self,'atoms') or self.atoms != atoms:
            self.atoms = atoms.copy()
            self.prepare_calculation()
            self.calculate()
            
    def prepare_calculation(self):
        raise NotImplementedException()

    def calculate(self):        
        pbc = self.atoms.get_pbc()
        if all(pbc):
            cell = self.atoms.get_cell()
        elif not any(pbc):
            # large enough cell for non-periodic calculation -
            # LAMMPS shrink-wraps automatically via input command 
            #       "periodic s s s"
            # below
            cell = 2 * np.max(np.abs(self.atoms.get_positions())) * np.eye(3)
        else: 
            print "WARNING: semi-periodic ASE cell detected -"
            print "         translation to proper LAMMPS input cell might fail"
            cell = self.atoms.get_cell()
        self.prism = prism(cell)
        self.run()

    def _lmp_alive(self):
        # Return True if this calculator is currently handling a running lammps process
        return self._lmp_handle and not isinstance(self._lmp_handle.poll(), int)

    def _lmp_end(self):
        # Close lammps input and wait for lammps to end. Return process return value
        if self._lmp_alive():
            self._lmp_handle.stdin.close()
            return self._lmp_handle.wait()
        
    def run(self):
        """Method which explicitely runs LAMMPS."""

        self.calls += 1

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


        # change into subdirectory for LAMMPS calculations
        cwd = os.getcwd()
        os.chdir(self.tmp_dir)
 

        # setup file names for LAMMPS calculation
        label = '%s%06d' % (self.label, self.calls)
        lammps_in = uns_mktemp(prefix='in_'+label, dir=self.tmp_dir)
        lammps_log = uns_mktemp(prefix='log_'+label, dir=self.tmp_dir)
        lammps_trj_fd = NamedTemporaryFile(prefix='trj_'+label, dir=self.tmp_dir,
                                           delete=(not self.keep_tmp_files))
        lammps_trj = lammps_trj_fd.name
        
        lammps_data_fd = NamedTemporaryFile(prefix='data_'+label, dir=self.tmp_dir,
                                            delete=(not self.keep_tmp_files))
        self.write_lammps_data(lammps_datafile=lammps_data_fd)
        lammps_data = lammps_data_fd.name
        lammps_data_fd.flush()


        # see to it that LAMMPS is started
        if not self._lmp_alive():
            # Attempt to (re)start lammps
            self._lmp_handle = Popen(lammps_cmd_line+lammps_options+['-log', '/dev/stdout'], 
                                    stdin=PIPE, stdout=PIPE)
        lmp_handle = self._lmp_handle


        # Create thread reading lammps stdout (for reference, if requested,
        # also create lammps_log, although it is never used)
        if self.keep_tmp_files:
            lammps_log_fd = open(lammps_log, 'w')
            fd = special_tee(lmp_handle.stdout, lammps_log_fd)
        else:
            fd = lmp_handle.stdout
        thr_read_log = Thread(target=self.read_lammps_log, args=(fd,))
        thr_read_log.start()


        # write LAMMPS input (for reference, also create the file lammps_in, 
        # although it is never used)
        if self.keep_tmp_files:
            lammps_in_fd = open(lammps_in, 'w')
            fd = special_tee(lmp_handle.stdin, lammps_in_fd)
        else:
            fd = lmp_handle.stdin
        self.write_lammps_in(lammps_in=fd, lammps_trj=lammps_trj, lammps_data=lammps_data)

        if self.keep_tmp_files:
            lammps_in_fd.close()


        # Wait for log output to be read (i.e., for LAMMPS to finish)
        # and close the log file if there is one
        thr_read_log.join()
        if self.keep_tmp_files:
            lammps_log_fd.close()

        if not self.keep_alive:
            self._lmp_end()

        exitcode = lmp_handle.poll()
        if exitcode and exitcode != 0:
            cwd = os.getcwd()
            raise RuntimeError('LAMMPS exited in %s with exit code: %d.' %\
                                   (cwd,exitcode))

        # A few sanity checks
        if len(self.thermo_content) == 0:
            raise RuntimeError('Failed to retreive any thermo_style-output')
        if int(self.thermo_content[-1]['atoms']) != len(self.atoms):
            # This obviously shouldn't happen, but if prism.fold_...() fails, it could
            raise RuntimeError('Atoms have gone missing')


        self.read_lammps_trj(lammps_trj=lammps_trj)
        lammps_trj_fd.close()
        lammps_data_fd.close()

        os.chdir(cwd)

        
    def write_lammps_data(self, lammps_datafile=None):
        """Method which writes a LAMMPS data file with atomic structure."""
        
        if (lammps_datafile == None):
            lammps_datafile = 'data.' + self.label
                          
        if isinstance(lammps_datafile, str):
            f = paropen(lammps_datafile, 'w')
            close_file = True
        else:
            # Presume fileobj acts like a fileobj
            f = lammps_datafile
            close_file = False
            
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

        if data.masses:
            f.write('Masses\n\n')
            for index, mass in enumerate(data.masses, 1):
                f.write('%i %f\n' % (index, mass))
            f.write('\n\n')
        
        def print_table(title, table):
            if len(table) == 0: return
            f.write('%s \n\n' % title)
            for index, row in enumerate(table, 1):
                f.write(('%d'+' %s'*len(row) +'\n') % ((index,) + tuple(row)))
            f.write('\n\n')
        
        for coeff_table in self.data.coeff_tables:
            print_table(coeff_table['title'], coeff_table['data'])
        
        # Atoms
        charges = self.atoms.get_charges()
        positions = p.pos_to_lammps(self.atoms.get_positions())
        if self.parameters.atom_style == 'full':
            table = [[tp, '1', a.charge, a.x, a.y, a.z]
                for tp, a in zip(self.data.atom_types, self.atoms)]
        elif  self.parameters.atom_style == 'charge':
            table = [[tp, a.charge, a.x, a.y, a.z]
                for tp, a in zip(self.data.atom_types, self.atoms)]
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
        
        if close_file:
            f.close()

    def write_lammps_in(self, lammps_in=None, lammps_trj=None, lammps_data=None):
        """Method which writes a LAMMPS in file with run parameters and settings."""

        if isinstance(lammps_in, str):
            f = paropen(lammps_in, 'w')
            close_in_file = True
        else:
            # Expect lammps_in to be a file-like object
            f = lammps_in
            close_in_file = False
            
        if self.keep_tmp_files:
            f.write('# (written by ASE)\n')

        # Write variables
        f.write('clear\n' +
                ('variable dump_file string "%s"\n' % lammps_trj) +
                ('variable data_file string "%s"\n' % lammps_data))

        parameters = self.parameters
        
        f.write('atom_style %s \n' % parameters.atom_style)
        
        pbc = self.atoms.get_pbc()
        f.write('units %s \n' % parameters.units)
        f.write('boundary %c %c %c \n' % tuple('sp'[x] for x in pbc))
        f.write('atom_modify sort 0 0.0 \n')
        if parameters.neighbor:
            f.write('neighbor %s \n' % (parameters.neighbor))
        if parameters.newton:
            f.write('newton %s \n' % (parameters.newton))
        f.write('\n')

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
        
        
        f.write('\n### read data \n')
        f.write('read_data %s\n' % lammps_data)
        
        # Extra pair coeffs
        for line in parameters.pair_coeffs:
            f.write('pair_coeff %s \n' % line)
   
        f.write('\n### run\n' +
                'fix fix_nve all nve\n' +
                ('dump dump_all all custom 1 %s id type x y z vx vy vz fx fy fz\n' % lammps_trj) )
        f.write(('thermo_style custom %s\n' +
                'thermo_modify flush yes\n' +
                'thermo 1\n') % (' '.join(self._custom_thermo_args)))

        if parameters.minimize:
            f.write('minimize %s\n' % parameters.minimize)
        if parameters.run:
            f.write('run %s\n' % parameters.run)
        if not (parameters.minimize or parameters.run):
            f.write('run 0\n')

        f.write('print "%s"\n' % CALCULATION_END_MARK)
        f.write('log /dev/stdout\n') # Force LAMMPS to flush log

        if close_in_file:
            f.close()

    def read_lammps_log(self, lammps_log=None, PotEng_first=False):
        """Method which reads a LAMMPS output log file."""

        if (lammps_log == None):
            lammps_log = self.label + '.log'

        if isinstance(lammps_log, str):
            f = paropen(lammps_log, 'w')
            close_log_file = True
        else:
            # Expect lammps_in to be a file-like object
            f = lammps_log
            close_log_file = False

        thermo_content = []
        line = f.readline()
        while line and line.strip() != CALCULATION_END_MARK:
            # get thermo output
            if line.startswith(self._custom_thermo_mark):
                m = True
                while m:
                    line = f.readline()
                    m = self._custom_thermo_re.match(line)
                    if m:
                        # create a dictionary between each of the thermo_style args
                        # and it's corresponding value
                        thermo_content.append(dict(zip(self._custom_thermo_args, 
                                                       map(float, m.groups()))))
            else:
                line = f.readline()

        if close_log_file:
            f.close()

        self.thermo_content = thermo_content

    def read_lammps_trj(self, lammps_trj=None, set_atoms=False):
        """Method which reads a LAMMPS dump file."""
        if (lammps_trj == None):
            lammps_trj = self.label + '.lammpstrj'

        f = paropen(lammps_trj, 'r')
        while True:
            line = f.readline()

            if not line:
                break

            #TODO: extend to proper dealing with multiple steps in one trajectory file
            if 'ITEM: TIMESTEP' in line:
                n_atoms = 0
                lo = [] ; hi = [] ; tilt = []
                id = [] ; type = []
                positions = [] ; velocities = [] ; forces = []

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
                # create corresponding index dictionary before iterating over atoms to (hopefully) speed up lookups...
                atom_attributes = {}
                for (i, x) in enumerate(line.split()[2:]):
                    atom_attributes[x] = i
                for n in range(n_atoms):
                    line = f.readline()
                    fields = line.split()
                    id.append( int(fields[atom_attributes['id']]) )
                    type.append( int(fields[atom_attributes['type']]) )
                    positions.append( [ float(fields[atom_attributes[x]]) for x in ['x', 'y', 'z'] ] )
                    velocities.append( [ float(fields[atom_attributes[x]]) for x in ['vx', 'vy', 'vz'] ] )
                    forces.append( [ float(fields[atom_attributes[x]]) for x in ['fx', 'fy', 'fz'] ] )
        f.close()

        # determine cell tilt (triclinic case!)
        if (len(tilt) >= 3):
            # for >=lammps-7Jul09 use labels behind "ITEM: BOX BOUNDS" to assign tilt (vector) elements ...
            if (len(tilt_items) >= 3):
                xy = tilt[tilt_items.index('xy')]
                xz = tilt[tilt_items.index('xz')]
                yz = tilt[tilt_items.index('yz')]
            # ... otherwise assume default order in 3rd column (if the latter was present)
            else:
                xy = tilt[0]
                xz = tilt[1]
                yz = tilt[2]
        else:
            xy = xz = yz = 0
        xhilo = (hi[0] - lo[0]) - xy - xz
        yhilo = (hi[1] - lo[1]) - yz
        zhilo = (hi[2] - lo[2])
    
        cell = [[xhilo,0,0],[xy,yhilo,0],[xz,yz,zhilo]]

        # assume that LAMMPS does not reorder atoms internally
        cell_atoms = np.array(cell)
        type_atoms = np.array(type)

        if self.atoms:
            cell_atoms = self.atoms.get_cell()

            # BEWARE: reconstructing the rotation from the LAMMPS output trajectory file
            #         fails in case of shrink wrapping for a non-periodic direction
            # -> hence rather obtain rotation from prism object used to generate the LAMMPS input
            #rotation_lammps2ase = np.dot(np.linalg.inv(np.array(cell)), cell_atoms)
            rotation_lammps2ase = np.linalg.inv(self.prism.R)

            type_atoms = self.atoms.get_atomic_numbers()
            positions_atoms = np.array( [np.dot(np.array(r), rotation_lammps2ase) for r in positions] )
            velocities_atoms = np.array( [np.dot(np.array(v), rotation_lammps2ase) for v in velocities] )
            forces_atoms = np.array( [np.dot(np.array(f), rotation_lammps2ase) for f in forces] )

        if (set_atoms):
            # assume periodic boundary conditions here (like also below in write_lammps)
            self.atoms = Atoms(type_atoms, positions=positions_atoms, cell=cell_atoms)

        self.forces = forces_atoms



class special_tee:
    """A special purpose, with limited applicability, tee-like thing.

    A subset of stuff read from, or written to, orig_fd, 
    is also written to out_fd.
    It is used by the lammps calculator for creating file-logs of stuff read from,
    or written to, stdin and stdout, respectively.
    """
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
        
    def is_skewed(self):
        acc = self.acc
        prism = self.get_lammps_prism()
        axy, axz, ayz = [np.abs(x) for x in prism[3:]]
        return (axy >= acc) or (axz >= acc) or (ayz >= acc)
        