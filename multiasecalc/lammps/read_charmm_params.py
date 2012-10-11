from lammpsbase import SequenceType, ImproperType, FFData

def read(infile):
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
		
		type = SequenceType(fields[0:2])
		parameters = [float(f) for f in fields[2:]]
		result.bond[type] = {'Bond Coeffs': parameters}
	
	# Angles
	for fields in parser(infile):
		if len(fields) == 5: fields += [0, 0]
		if len(fields) != 7: continue
		
		type = SequenceType(fields[0:3])
		parameters = [float(f) for f in fields[3:]]
		result.angle[type] = {'Angle Coeffs': parameters}
		
	# Dihedrals
	for fields in parser(infile):
		if len(fields) == 7: fields += [0]
		if len(fields) != 8: continue
		
		type = SequenceType(fields[0:4])
		parameters = [float(fields[4]), int(fields[5]), int(float(fields[6])), float(fields[7])]
		result.dihedral[type] = {'Dihedral Coeffs': parameters}
		
	# Impropers
	for fields in parser(infile):
		if len(fields) != 7: continue
		
		type = ImproperType(central_type=fields[0], other_types=fields[1:4])
		parameters = [float(f) for f in fields[4:7:2]]
		result.improper[type] = {'Improper Coeffs': parameters}
	
	# Nonbonded
	for fields in parser(infile):
		if len(fields) == 4: fields += [0, 0, 0]
		if len(fields) != 7: continue
		
		eps = -float(fields[2])
		rmin = float(fields[3])*2
		sigma = rmin*2**(-1./6)
		
		type = fields[0]
		parameters = [eps, sigma]
		result.atom[type] = {'Pair Coeffs': parameters}
		
	return result
	