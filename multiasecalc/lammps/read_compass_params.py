from ffdata import SequenceType, ImproperType, FFData
import numpy as np
import re

def fileIterator(infile):
	data  = infile.read()
	data = re.sub(r'\s*![^\n]*\n', '\n', data)
	
	frcTableExp = re.compile(r"""
		\#(\S+)\s+\w+\s*\n+
		(?:>[^\n]*\n)*
		\s*
		(?:@[^\n]*\n)*
		\s*
		(.+?)\n\n
		""", re.VERBOSE | re.DOTALL)

	for match in frcTableExp.finditer(data):
		yield match.groups()



# ******************************************

def read(infile):
	ff_data = FFData()
	ff_data.class2 = True
	
	R0 = {}
	theta0 = {}
	
	def get_R0(pair):
		for tp in SequenceType(pair).variations():
			R = R0.get(tp)
			if R != None: return R
		#print 'No R0 for %s' % (pair,)

	def get_theta0(triplet):
		for tp in SequenceType(triplet).variations():
			th = theta0.get(tp)
			if th != None: return th
		#print 'No theta0 for %s' % (triplet,)
	
	it = fileIterator(infile)
	
	# This is a hack... In pcff.frc the wilson_out_of_plane table comes twice, and
	# on the second time it's for cff91_auto. We want to ignore that.
	#tables = dict(it)
	tables = {}
	for title, table in it:
		if title not in tables:
			tables[title] = table
	
	# **** Equivalence ****
	table = tables['equivalence']
	for row in table.splitlines():
		fields = row.split()[2:]
		order = ['atom', 'bond', 'angle', 'dihedral', 'improper']
		ff_data.add_equivalence(fields[0], **dict(zip(order, fields[1:])))
	
	def iterateTable(data, category):
		nAtomTypes = dict(NonB=1, Bond=2, Angle=3, Torsion=4, OOP=4)[category]
		for line in data.splitlines():
			row = line.split()
			actual_types = tuple(row[2:2+nAtomTypes])
			values = [float(val) for val in row[2+nAtomTypes:]]
			yield actual_types, values
	
	# **** Bond coeffs ****
	table = tables['quartic_bond']
	for atomTypes, values in iterateTable(table, 'Bond'):
		type = SequenceType(atomTypes)
		ff_data.add('bond', type, 'Bond Coeffs', values)
		R0[atomTypes] = values[0]
		try:
			dihedral_type = ff_data.get_actual_type('dihedral', type)
			if dihedral_type not in R0:
				R0[dihedral_type] = values[0]
		except KeyError: pass
	
	# **** Angle coeffs ****
	table = tables['quartic_angle']
	for atomTypes, values in iterateTable(table, 'Angle'):
		type = SequenceType(atomTypes)
		ff_data.add('angle', type, 'Angle Coeffs', values)
		theta0[type] = values[0]
		try:
			improper_type = ff_data.get_actual_type('improper', type)
			if improper_type not in theta0:
				theta0[improper_type] = values[0]
		except KeyError: pass
	 
	# **** Bond/bond ****
	table = tables['bond-bond']
	for types, values in iterateTable(table, 'Angle'):
		R1 = get_R0(types[:2])
		R2 = get_R0(types[1:])
		if None in (R1, R2): continue
		ff_data.add('angle', SequenceType(types), 'BondBond Coeffs', [values[0], R1, R2])
	
	# **** Bond/angle ****
	table = tables['bond-angle']
	for types, values in iterateTable(table, 'Angle'):
		R1 = get_R0(types[:2])
		R2 = get_R0(types[1:])
		if None in (R1, R2): continue
		K1 = values[0]
		K2 = K1
		if len(values) == 2: K2 = values[1]
		type = SequenceType(types)
		ff_data.add('angle', type, 'BondAngle Coeffs', [ K1, K2, R1, R2])
		if not 'BondBond Coeffs' in ff_data.get_params('angle', type):
			# In practice LAMMPS uses the R1 and R2 values from BondBond coeffs so we
			# need to make sure those are given
			ff_data.add('angle', type, 'BondBond Coeffs', [0.0, R1, R2])
	
	# **** Bond/bond 1-3 ***
	table = tables['bond-bond_1_3']
	for types, values in iterateTable(table, 'Torsion'):
		R1 = get_R0(types[:2])
		R2 = get_R0(types[2:])
		if None in (R1, R2): continue
		ff_data.add('dihedral', SequenceType(types), 'BondBond13 Coeffs', [values[0], R1, R2])

	# **** Torsion coeffs ****
	table = tables['torsion_3']
	for types, values in iterateTable(table, 'Torsion'):
		ff_data.add('dihedral', SequenceType(types, wildcard='*'), 'Dihedral Coeffs', values)

	# **** End bond/torsion ***
	table = tables['end_bond-torsion_3']
	for types, values in iterateTable(table, 'Torsion'):
		R1 = get_R0(types[:2])
		R2 = get_R0(types[2:])
		if None in (R1, R2): continue
		B = values[:3]
		C = B
		if len(values) > 3: C = values[3:]
		ff_data.add('dihedral', SequenceType(types), 'EndBondTorsion Coeffs',  B + C + [R1, R2])
	
	# **** Middle bond/torsion ****
	table = tables['middle_bond-torsion_3']
	for types, values in iterateTable(table, 'Torsion'):
		R2 = get_R0(types[1:3])
		if not R2: continue
		ff_data.add('dihedral', SequenceType(types), 'MiddleBondTorsion Coeffs', values + [R2])
		
	# **** Angle/torsion ****
	table = tables['angle-torsion_3']
	for types, values in iterateTable(table, 'Torsion'):
		theta1 = get_theta0(types[:3])
		theta2 = get_theta0(types[1:])
		if None in (theta1, theta2): continue
		D  = values[:3]
		E = D
		if len(values) > 3: E = values[3:]
		ff_data.add('dihedral', SequenceType(types), 'AngleTorsion Coeffs', D + E + [theta1, theta2])
		
	
	# **** Wilson out of plane ****
	table = tables['wilson_out_of_plane']
	for types, values in iterateTable(table, 'OOP'):
		ff_data.add('improper', ImproperType(types, class2=True), 'Improper Coeffs', values)

		
	# ****** angle/angle ******
	table = tables['angle-angle']
	aaCoeffs = {}
	aa_dihedrals = set()
	for types, values in iterateTable(table, 'OOP'):
		aaCoeffs[ (types[1], types[2], (types[0], types[3])) ] = values[0]
		aa_dihedrals.add(ImproperType(atom_types=types, class2=True))

	for type in aa_dihedrals:
		ind, actual_type = ff_data.find('improper', range(4), type)
		if actual_type:
			type = actual_type
		
		j, (i, k, l) = type.get_types()
		
		key1 = (j, k, (i, l))
		if not key1 in aaCoeffs: key1 = (j, k, (l, i))
		M1 = aaCoeffs.get(key1)
		
		key2 = (j, i, (k, l))
		if not key2 in aaCoeffs: key2 = (j, i, (l, k))
		M2 = aaCoeffs.get(key2)
		
		key3 = (j, l, (k, i))
		if not key3 in aaCoeffs: key3 = (j, l, (i, k))
		M3 = aaCoeffs.get(key3)
		
		if None in (M1, M2, M3):
			#print 'Angle-angle: not all combinations found for', i, j, k, l
			continue
		
		th1 = get_theta0([i,j,k])
		th2 = get_theta0([i,j,l])
		th3 = get_theta0([k,j,l])
		if None in (th1, th2, th3): continue
		ff_data.add('improper', type, 'AngleAngle Coeffs', [M1, M2, M3, th1, th2, th3])
		

	# ****** angle/angle/torsion ******
	table = tables['angle-angle-torsion_1']
	for types, values in iterateTable(table, 'Torsion'):
		theta1 = get_theta0(types[:3])
		theta2 = get_theta0(types[1:])
		if None in (theta1, theta2): continue
		v = [values[0], theta1, theta2]
		ff_data.add('dihedral', SequenceType(types), 'AngleAngleTorsion Coeffs', v)

	# ****** nonbonded ******
	table = tables['nonbond(9-6)']
	for types, values in iterateTable(table, 'NonB'):
		ff_data.add('atom', types[0], 'Pair Coeffs', values[::-1])
		
	# ****** bond increments ******
	bond_increments = {}
	table = tables['bond_increments']
	for types, values in iterateTable(table, 'Bond'):
		bond_increments[types] = values
		
	return ff_data, bond_increments
