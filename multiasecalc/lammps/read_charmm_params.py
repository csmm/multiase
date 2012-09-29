def read(infile):
	class FFData: 
		def __init__(self):
			#self.mass = {}	
			self.bond = {}
			self.angle = {}
			self.dihedral = {}
			self.improper = {'empty': [0, 0]}
			self.nonbonded = {}

	result = FFData()
	
	def parser(infile):
		for line in infile:
			if line.startswith('BONDS') or \
				line.startswith('ANGLES') or \
				line.startswith('DIHEDRALS') or \
				line.startswith('IMPROPERS') or \
				line.startswith('NONBONDED') or \
				line.startswith('HBOND'):
				break
			
			line = line.replace('!HGTIP', 'HGTIP')
			line = line.replace('!OGTIP', 'OGTIP')
			yield line.split('!')[0].split()
	
	# Skip masses
	for item in parser(infile): pass
		
	# Bonds
	for fields in parser(infile):
		if len(fields) != 4: continue
		result.bond[tuple(fields[0:2])] = [float(f) for f in fields[2:]]
	
	# Angles
	for fields in parser(infile):
		if len(fields) == 5: fields += [0, 0]
		if len(fields) != 7: continue
		
		result.angle[tuple(fields[0:3])] = [float(f) for f in fields[3:]]
		
	# Dihedrals
	for fields in parser(infile):
		if len(fields) == 7: fields += [0]
		if len(fields) != 8: continue
		
		result.dihedral[tuple(fields[0:4])] = [float(fields[4]), int(fields[5]), int(float(fields[6])), float(fields[7])]
		
	# Impropers
	for fields in parser(infile):
		if len(fields) != 7: continue
		result.improper[tuple(fields[0:4])] = [float(f) for f in fields[4:7:2]]
	
	# Nonbonded
	for fields in parser(infile):
		if len(fields) == 4: fields += [0, 0, 0]
		if len(fields) != 7: continue
		
		eps = -float(fields[2])
		rmin = float(fields[3])*2
		sigma = rmin*2**(-1./6)
		
		result.nonbonded[(fields[0],)] = [eps, sigma]
		
	return result
	
if __name__ == '__main__':
	result = read(open('par_all36_cgenff.prm'))
	print result.nonbonded