
funcDict = {}

def getType(atom, neighbors):
	return funcDict[atom.symbol]([n.symbol for n in neighbors])

# **** Carbon ****
def C(neighbors):
	if len(neighbors) == 3:
		if neighbors.count('C') == 2:
			return 'c3a'
			
		if 'O' in neighbors:
			return "c3'"
	
	if neighbors.count('H') >= 2:
		if 'O' in neighbors:
			return 'c4o'
		else: return 'c4'
	
	if neighbors.count('H') == 1:
		return 'c43'
	
	# 4 heavy atoms
	return 'c44'
funcDict['C'] = C

# **** Hydrogen ****
def H(neighbors):
	if 'H' in neighbors:
		return 'h1h'
	
	# default
	return 'h1'
funcDict['H'] = H
	
# **** Oxygen ****
def O(neighbors):
	#if list(neighbors) == ['H', 'H']:
	#	return 'o2*'
		
	if list(neighbors) == ['C', 'C']:
		return 'o2e'
	
	if neighbors.count('H') == 1:
		return 'o2h'
		
	if len(neighbors) == 1:
		return 'o1='
	
	# default
	return 'o2'
funcDict['O'] = O

# **** Nitrogen ****
def N(neighbors):
	if list(neighbors) == ['O', 'O']:
		return 'n2o'
		
	if list(neighbors) == ['N']:
		return 'n1n'
		
	if list(neighbors) == ['O']:
		return 'n1o'
		
	# default
	return 'n2='
funcDict['N'] = N

'''
class eqTableEntry:
	def __init__(self, nonb, bond, angl, tors, oop):
		self.nonbonded=nonb; self.bonded=bond; self.angle=angl; self.torsion=tors; self.outOfPlane=oop

equivalenceTable = {
#  formalType (nonb. bond  angl  tors. oop)
	'c4' : eqTableEntry('c4', 'c4', 'c4', 'c4', 'c4'),
	'c43': eqTableEntry('c43', 'c4', 'c4', 'c4', 'c4'),
	'c44': eqTableEntry('c44', 'c4', 'c4', 'c4', 'c4'),
	'c4o': eqTableEntry('c4o', 'c4', 'c4', 'c4', 'c4'),
	'h1' : eqTableEntry('h1', 'h1', 'h1', 'h1', 'h1'),
	'h1h': eqTableEntry('h1h', 'h1h', 'h1', 'h1', 'h1'),
	'o2' : eqTableEntry('o2', 'o2', 'o2', 'o2', 'o2'),
	'o2h' : eqTableEntry('o2h', 'o2h', 'o2', 'o2', 'o2'),
	'o2*': eqTableEntry('o2*', 'o2*', 'o2*', 'o2*', 'o2*'),
}
'''
