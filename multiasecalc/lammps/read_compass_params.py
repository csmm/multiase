from lammpsbase import SequenceType, ImproperType, FFData
import numpy as np

# itertools.product not in python 2.4

#from itertools import product

def product(*args, **kwds):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = map(tuple, args) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)



def fileIterator(infile):
	import re
	frcTableExp = re.compile(r"""
		\#(\S+)\s+\w+\s*\n+
		(?:>[^\n]*\n)*
		\s*
		(?:@[^\n]*\n)*
		\s*
		(?:![^\n]*\n)*
		(.+?)\n\n
		""", re.VERBOSE | re.DOTALL)

	for match in frcTableExp.finditer(infile.read()):
		yield match.groups()



# ******************************************

def read(infile):
	ff_data = FFData()
	ff_data.class2 = True
	
	R0 = {}
	theta0 = {}

	it = fileIterator(infile)
	
	tables = dict(it)
	
	# **** Equivalence ****
	table = tables['equivalence']
	equivalence = np.array([row.split()[2:] for row in table.split('\n')])
	
	def formalTypes(actualTypes, category):
		columns = ['NonB', 'Bond', 'Angle', 'Torsion', 'OOP']
		col = columns.index(category)+1
		def formalTypesSingle(actualType):
			indices = np.where(equivalence[:,col] == actualType)[0]
			return equivalence[indices, 0]
			
		for types in product(*(formalTypesSingle(at) for at in actualTypes)):
			yield types
	
	def iterateTable(data, category):
		nAtomTypes = dict(NonB=1, Bond=2, Angle=3, Torsion=4, OOP=4)[category]
		for line in data.split('\n'):
			if '*' in line: continue # Ignore these for now
			row = line.split()
			actualTypes = row[2:2+nAtomTypes]
			values = [float(val) for val in row[2+nAtomTypes:]]
			for types in formalTypes(actualTypes, category):
				yield types, values
	
	# **** Bond coeffs ****
	table = tables['quartic_bond']
	for atomTypes, values in iterateTable(table, 'Bond'):
		type = SequenceType(atomTypes)
		ff_data.bond[type] = {'Bond Coeffs': values}
		R0[type] = values[0]
	
	# **** Angle coeffs ****
	table = tables['quartic_angle']
	ff_data.bond['Angle Coeffs'] = {}
	for atomTypes, values in iterateTable(table, 'Angle'):
		type = SequenceType(atomTypes)
		ff_data.angle[type] = {'Angle Coeffs': values}
		theta0[type] = values[0]
	 
	# **** Bond/bond ****
	table = tables['bond-bond']
	for types, values in iterateTable(table, 'Angle'):
		bond1 = SequenceType(types[:2])
		bond2 = SequenceType(types[1:])
		if not bond1 in R0 or not bond2 in R0: continue
		ff_data.add('angle', SequenceType(types), 'BondBond Coeffs', [values[0], R0[bond1], R0[bond2]])
	
	# **** Bond/bond 1-3 ***
	table = tables['bond-bond_1_3']
	for types, values in iterateTable(table, 'Torsion'):
		bond1 = SequenceType(types[:2])
		bond2 = SequenceType(types[2:])
		if not bond1 in R0 or not bond2 in R0: continue
		ff_data.add('dihedral', SequenceType(types), 'BondBond13 Coeffs', [values[0], R0[bond1], R0[bond2]])
		
	# **** Bond/angle ****
	table = tables['bond-angle']
	for types, values in iterateTable(table, 'Angle'):
		bond1 = SequenceType(types[:2])
		bond2 = SequenceType(types[1:])
		if not bond1 in R0 or not bond2 in R0: continue
		K2 = 0
		if len(values) == 2: K2 = values[1]
		ff_data.add('angle', SequenceType(types), 'BondAngle Coeffs', [ values[0], K2, R0[bond1], R0[bond2]])

	# **** Torsion coeffs ****
	table = tables['torsion_3']
	for types, values in iterateTable(table, 'Torsion'):
		ff_data.add('dihedral', SequenceType(types), 'Dihedral Coeffs', values)

	# **** End bond/torsion ***
	table = tables['end_bond-torsion_3']
	for types, values in iterateTable(table, 'Torsion'):
		bond1 = SequenceType(types[:2])
		bond2 = SequenceType(types[2:])
		if not bond1 in R0 or not bond2 in R0: continue
		B = values[:3]
		C = B
		if len(values) > 3: C = values[3:]
		ff_data.add('dihedral', SequenceType(types), 'EndBondTorsion Coeffs',  B + C + [R0[bond1], R0[bond2]])
	
	# **** Middle bond/torsion ****
	table = tables['middle_bond-torsion_3']
	for types, values in iterateTable(table, 'Torsion'):
		bond2 = SequenceType(types[1:3])
		if not bond2 in R0: continue
		ff_data.add('dihedral', SequenceType(types), 'MiddleBondTorsion Coeffs', values + [R0[bond2]])
		
	# **** Angle/torsion ****
	table = tables['angle-torsion_3']
	for types, values in iterateTable(table, 'Torsion'):
		angle1 = SequenceType(types[:3])
		angle2 = SequenceType(types[1:])
		if not angle1 in theta0 or not angle2 in theta0: continue
		D  = values[:3]
		E = D
		if len(values) > 3: E = values[3:]
		ff_data.add('dihedral', SequenceType(types), 'AngleTorsion Coeffs', D + E + [theta0[angle1], theta0[angle2]])
		
	
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
		aa_dihedrals.add(ImproperType(types, class2=True))

	for type in aa_dihedrals:
		j, (i, k, l) = type.get_types()
		
		key1 = (j, k, (i, l))
		if not key1 in aaCoeffs:
			key1 = (j, i, (l, i))
			if not key1 in aaCoeffs: continue
		M1 = aaCoeffs.get(key1)
		
		key2 = (j, i, (k, l))
		if not key2 in aaCoeffs: key2 = (j, i, (l, k))
		M2 = aaCoeffs.get(key2)
		
		key3 = (j, l, (k, i))
		if not key3 in aaCoeffs: key3 = (j, l, (i, k))
		M3 = aaCoeffs.get(key3)
		
		if None in (M1, M2, M3):
			print 'Angle-angle: not all combinations found for', i, j, k, l
			continue
		
		th1 = theta0.get(SequenceType([i,j,k]))
		th2 = theta0.get(SequenceType([i,j,l]))
		th3 = theta0.get(SequenceType([k,j,l]))
		if None in (th1, th2, th3): continue
		ff_data.add('improper', type, 'AngleAngle Coeffs', [M1, M2, M3, th1, th2, th3])
		

	# ****** angle/angle/torsion ******
	table = tables['angle-angle-torsion_1']
	for types, values in iterateTable(table, 'Torsion'):
		angle1 = SequenceType(types[:3])
		angle2 = SequenceType(types[1:])
		if not angle1 in theta0 or not angle2 in theta0: continue
		v = [values[0], theta0[angle1], theta0[angle2]]
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
