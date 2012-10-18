from ase.atoms import Atoms
from ase.quaternions import Quaternions
from ase.parallel import paropen
import numpy as np

""" This could be later merged with ase.io.lammps """

def read_lammps_dump(fileobj, index=-1):
	"""Method which reads a LAMMPS dump file."""
	if isinstance(fileobj, str):
		f = paropen(fileobj)
	else:
		f = fileobj

	# load everything into memory
	lines = f.readlines()

	natoms = 0
	images = []

	while len(lines) > natoms:
		line = lines.pop(0)

		if 'ITEM: TIMESTEP' in line:
			natoms = 0
			lo = [] ; hi = [] ; tilt = []
			#id = [] ; types = []
			positions = None
			velocities = None
			forces = None
			quaternions = None
			charges = None

		if 'ITEM: NUMBER OF ATOMS' in line:
			line = lines.pop(0)
			natoms = int(line.split()[0])
			
		if 'ITEM: BOX BOUNDS' in line:
			# save labels behind "ITEM: BOX BOUNDS" in triclinic case (>=lammps-7Jul09)
			tilt_items = line.split()[3:]
			for i in range(3):
				line = lines.pop(0)
				fields = line.split()
				lo.append(float(fields[0]))
				hi.append(float(fields[1]))
				if (len(fields) >= 3):
					tilt.append(float(fields[2]))

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
			xhilo = (hi[0] - lo[0]) - abs(xy) - abs(xz)
			yhilo = (hi[1] - lo[1]) - abs(yz)
			zhilo = (hi[2] - lo[2])
			celldispx = lo[0] - min(0, xy) - min(0, xz)
			celldispy = lo[1] - min(0, yz)
			celldispz = lo[2]

			cell = [[xhilo,0,0],[xy,yhilo,0],[xz,yz,zhilo]]
			celldisp = [[celldispx, celldispy, celldispz]]
				
		if 'ITEM: ATOMS' in line:
			atom_attributes = {}
			for (i, x) in enumerate(line.split()[2:]):
				atom_attributes[x] = i
			
			positions = np.zeros((natoms, 3))*np.nan
			types = np.empty(natoms, dtype=int)
			if 'vx' in atom_attributes:
				velocities = np.zeros((natoms, 3))*np.nan
			
			if 'fx' in atom_attributes:
				forces = np.zeros((natoms, 3))*np.nan
			
			if 'c_q[1]' in atom_attributes:
				quaternions = np.zeros((natoms, 4))*np.nan
			
			if 'q' in atom_attributes:
				charges = np.zeros(natoms)
			
			def get_values(fields, keys):
				return [float(fields[atom_attributes[key]]) for key in keys]
			
			for n in range(natoms):
				line = lines.pop(0)
				fields = line.split()
				id = int(fields[atom_attributes['id']])
				#id.append( int(fields[atom_attributes['id']]) )
				types[id-1] = int(fields[atom_attributes['type']])
				try:
					positions[id-1] = get_values(fields, ['x', 'y', 'z'])
				except:
					print id, positions.shape, n_atoms, line
					raise
				if velocities != None:
					velocities[id-1] = get_values(fields, ['vx', 'vy', 'vz'])
				if forces != None:
					forces[id-1] = get_values(fields, ['fx', 'fy', 'fz'])
				if quaternions != None:
					quaternions[id-1] = get_values(fields, ['c_q[1]', 'c_q[2]', 'c_q[3]', 'c_q[4]'])
				if charges != None:
					charges[id-1] = get_values(fields, 'q')[0]

			#if len(quaternions):
			#    images.append(Quaternions(symbols=types,
			#                              positions=positions,
			#                              cell=cell, celldisp=celldisp,
			#                              quaternions=quaternions, info=info))

			atoms = Atoms(symbols=types, positions=positions, cell=cell)
			if velocities != None: 
				atoms.set_velocities(velocities)
			if charges != None:
				atoms.set_charges(charges)
			info = dict(celldisp=celldisp)
			if forces != None:
				info['forces'] = forces
			atoms.info = info
			images.append(atoms)

	return images[index]

