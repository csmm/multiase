from itertools import product

def itemize(data, nAtomTypes):
	for row in data.split('\n'):
		row = row.split()
		yield tuple(row[2:2+nAtomTypes]), [float(val) for val in row[2+nAtomTypes:]]

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

	title = yield
	for match in frcTableExp.finditer(infile.read()):
		if match.group(1) == title:
			title = yield match.group(2)



# ******************************************

def read(infile):
	''' The main function'''
	
	class ffData: 
		def __init__(self):
			self.mass = {}
			self.equivalence = {}
			self.bond = {}
			self.angle = {}
			self.bondBond = dict(empty = [0]*3)
			self.bondBond13 = dict(empty = [0]*3)
			self.bondAngle = dict(empty = [0]*3)
			self.torsion = {}
			self.endBondTorsion = {}
			self.middleBondTorsion = {}
			self.angleTorsion = {}
			self.improper = dict(empty = [0]*2)
			self.angleAngle = dict(empty = [0]*4)
			self.angleAngleTorsion = {}
			self.pair = {}
			self.bondIncrement = {}
			
	result = ffData()
	
	R0 = {}
	theta0 = {}
	equivalence = result.equivalence

	it = fileIterator(infile)
	it.next()
	
	# **** Masses ****
	table = it.send('atom_types')
	for row in table.split('\n'):
		row = row.split()
		result.mass[row[2]] = float(row[3])
	
	# **** Equivalence ****
	table = it.send('equivalence')
	for atomTypes, values in itemize(table, 6):
		nb, b, a, t, oop = atomTypes[1:]
		result.equivalence[atomTypes[0]] = dict(nonbonded=nb, bonded=b, angle=a, torsion=t, outOfPlane=oop)
	
	def typeCombinations(int1_types, interaction1, interaction2):
		''' Returns an iterable for going through all combinations of possible interaction2 type
			combinations for the given set of interaction1 types. Example:
			interaction1 = 'angle'
			int1_types = ('c4', 'o2', 'c4')
			interaction2 = 'bonded'
			returns: ('c4', 'o2', 'c4'), ('c4', 'o2e', 'c4'), ('c43', 'o2h', 'c4') '''
		
		def possibleTypes(int1_type, interaction1, interaction2):
			return set(e[interaction2] for e in equivalence.values() if e[interaction1] == int1_type)
		
		pt = ( possibleTypes(t, interaction1, interaction2) for t in int1_types )
		for types in product(*pt):
			yield types
	
	# **** Bond coeffs ****
	table = it.send('quartic_bond')
	for atomTypes, values in itemize(table, 2):
		R0[atomTypes] = values[0]
		result.bond[atomTypes] = values

	def getR0(types):
		''' Helper function '''
		if types not in R0: types = types[::-1]
		return R0.get(types)
		
	# **** Angle coeffs ****
	table = it.send('quartic_angle')
	for atomTypes, values in itemize(table, 3):
		theta0[atomTypes] = values[0]
		result.angle[atomTypes] = values
		
	def getTheta0(types):
		''' Helper function '''
		if types not in theta0: types = types[::-1]
		return theta0.get(types)
			
	# **** Bond/bond ****
	table = it.send('bond-bond')
	for atomTypes, values in itemize(table, 3):
		for types in typeCombinations(atomTypes, 'angle', 'bonded'):
			r1 = getR0(types[:2])
			r2 = getR0(types[1:])
			if None not in (r1, r2):
				result.bondBond[types] = [ values[0], r1, r2 ]
	
	# **** Bond/bond 1-3 ***
	table = it.send('bond-bond_1_3')
	for atomTypes, values in itemize(table, 4):
		for types in typeCombinations(atomTypes, 'torsion', 'bonded'):
			r1 = getR0(types[:2])
			r2 = getR0(types[2:])
			if None not in (r1, r2):
				result.bondBond13[types] = [ values[0], r1, r2]
		
	# **** Bond/angle ****
	table = it.send('bond-angle')
	for atomTypes, values in itemize(table, 3):
		for types in typeCombinations(atomTypes, 'angle', 'bonded'):
			r1 = getR0(types[:2])
			r2 = getR0(types[1:])
			if None not in (r1, r2):
				K2 = 0
				if len(values) == 2: K2 = values[1]
				result.bondAngle[types] = [ values[0], K2, r1, r2 ]

	# **** Torsion coeffs ****
	table = it.send('torsion_3')
	for atomTypes, values in itemize(table, 4):
		result.torsion[atomTypes] = values

	# **** End bond/torsion ***
	table = it.send('end_bond-torsion_3')
	for atomTypes, values in itemize(table, 4):
		for types in typeCombinations(atomTypes, 'angle', 'bonded'):
			r1 = getR0(types[:2])
			r2 = getR0(types[2:])
			if None not in (r1, r2):
				B = values[:3]
				C = B
				if len(values) > 3: C = values[3:]
				result.endBondTorsion[types] = B + C + [r1, r2]
	
	# **** Middle bond/torsion ****
	table = it.send('middle_bond-torsion_3')
	for atomTypes, values in itemize(table, 4):
		for types in typeCombinations(atomTypes, 'angle', 'bonded'):
			r2 = getR0(types[1:3])
			if r2 != None:
				result.middleBondTorsion[types] = values + [r2]
		
	# **** Angle/torsion ****
	table = it.send('angle-torsion_3')
	for atomTypes, values in itemize(table, 4):
		for types in typeCombinations(atomTypes, 'torsion', 'angle'):
			th1 = getTheta0(types[:3])
			th2 = getTheta0(types[1:])
			if None not in (th1, th2):
				D  = values[:3]
				E = D
				if len(values) > 3: E = values[3:]
				result.angleTorsion[types] = D + E + [th1, th2]
	
	# **** Wilson out of plane ****
	table = it.send('wilson_out_of_plane')
	for atomTypes, values in itemize(table, 4):
		result.improper[atomTypes] = values

		
	# ****** angle/angle ******
	table = it.send('angle-angle')
	aaCoeffs = {}
	for atomTypes, values in itemize(table, 4):
		aaCoeffs[ (atomTypes[1], atomTypes[2], (atomTypes[0], atomTypes[3])) ] = values[0]

	alreadyDone = []
	for key, M1 in aaCoeffs.items():
		if key in alreadyDone: continue
		
		j, k, (i, l) = key
		key2 = (j, i, (k, l))
		if not key2 in aaCoeffs: key2 = (j, i, (l, k))
		M2 = aaCoeffs.get(key2)
		
		key3 = (j, l, (k, i))
		if not key3 in aaCoeffs: key3 = (j, l, (i, k))
		M3 = aaCoeffs.get(key3)
		
		alreadyDone += [key, key2, key3]
		if None in (M2, M3):
			print 'Angle-angle: not all combinations found for', i, j, k, l
			continue
		
		for a, b, c, d in typeCombinations((i, j, k, l), 'torsion', 'angle'):
			th1 = getTheta0((a, b, c))
			th2 = getTheta0((a, b, d))
			th3 = getTheta0((c, b, d))
			if None not in (th1, th2, th3):
				result.angleAngle[(a, b, c, d)] = [M1, M2, M3, th1, th2, th3]
		

	# ****** angle/angle/torsion ******
	table = it.send('angle-angle-torsion_1')
	for atomTypes, values in itemize(table, 4):
		for types in typeCombinations(atomTypes, 'torsion', 'angle'):
			th1 = getTheta0(atomTypes[:3])
			th2 = getTheta0(atomTypes[1:])
			if None not in (th1, th2):
				result.angleAngleTorsion[types] = [values[0], th1, th2]

	# ****** nonbonded ******
	table = it.send('nonbond(9-6)')
	for atomTypes, values in itemize(table, 1):
		result.pair[atomTypes] = values[::-1]
		
	# ****** bond increments ******
	table = it.send('bond_increments')
	for atomTypes, values in itemize(table, 2):
		result.bondIncrement[atomTypes] = values
		
	return result

	
'''
infile = open('/home/tuukka/Documents/compass.frc')

result = read(infile)

for key, value in result.items():
	print key
	for k, v in value.items():
		print k, ':', v
	print

#for e in angle_angle: print e
'''
