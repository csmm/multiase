from ffdata import FFData, SequenceType

nm = 10
kJmol = 0.239005736

def read(files, return_charges=False):
	ffdata = FFData()
	if return_charges:
		charges = {}
	else: charges = None
	for filename in files:
		f = open(filename)
		read_file(f, ffdata, charges)
		f.close()
		
	if return_charges: return ffdata, charges
	else: return ffdata
	
def read_file(infile, ffdata, charge_dict=None):
	state = None
	for line in infile:
		line = line.split(';')[0].strip()
		
		if not line:
			continue
		if line[0] == '[':
			state = line[1:-1].strip()
			continue
		if line[0] == '#':
			continue
		
		fields = line.split()
		
		if state == 'atomtypes':
			add_nonbonded(fields, ffdata, charge_dict)
		elif state == 'bondtypes':
			add_bond(fields, ffdata)
		elif state == 'angletypes':
			add_angle(fields, ffdata)
		elif state == 'dihedraltypes':
			add_dihedral(fields, ffdata)


#; name      at.num  mass     charge ptype  sigma      epsilon
#H0           1       1.008   0.0000  A   2.47135e-01  6.56888e-02

#; name  bond_type    mass    charge   ptype          sigma      epsilon
# opls_001   C	6      12.01100     0.500       A    3.75000e-01  4.39320e-01 ; SIG

def add_nonbonded(fields, ffdata, charge_dict):
	epsilon = float(fields[-1]) * kJmol
	sigma   = float(fields[-2]) * nm
	charge  = float(fields[-4])
	name = fields[0]
	if len(fields) == 8:
		bondtype = fields[1]
		ffdata.add_equivalence(name, atom=name, bond=bondtype, angle=bondtype, dihedral=bondtype)
	ffdata.add('atom', name, 'Pair Coeffs', [epsilon, sigma])
	if charge_dict != None:
		charge_dict[name] = charge


#; i    j  func       b0          kb
#  CT H0         1    0.10900   284512.0

def add_bond(fields, ffdata):
	type = SequenceType(fields[:2])
	K = 0.5*float(fields[-1]) * kJmol/nm**2
	r0 = float(fields[-2]) * nm
	ffdata.add('bond', type, 'Bond Coeffs', [K, r0])


#;  i    j    k  func       th0       cth
#H0  CT  H0           1   109.500    292.880

def add_angle(fields, ffdata):
	type = SequenceType(fields[:3])
	K = float(fields[-1]) * kJmol
	th0 = float(fields[-2])
	ffdata.add('angle', type, 'Angle Coeffs', [K, th0])


#;i  j   k  l	 func      phase      kd      pn
#CA  CA  CA  OH       4      180.00     4.60240     2 

#;  i    j    k    l   func     coefficients
#  Br     C      CB     CT      3      0.00000   0.00000   0.00000   0.00000   0.00000   0.00000

def add_dihedral(fields, ffdata):
	for name in fields[:4]:
		# Check for weird entries (I guess X name name X can be abbreviated with name name)
		if not name[0].isalpha():
			return
			
	type = SequenceType(fields[:4], wildcard='X')
	if fields[4] in ('4', '9'):
		K = float(fields[-2]) * kJmol
		n = int(fields[-1])
		d = int(float(fields[-3]))
		parameters = [K, n, d, 0]
	elif fields[4] == '3':
		# C5 should be zero
		if float(fields[10]) != 0:
			print 'RB dihedral has an unsupported nonzero C5 coeff. Ignoring...'
			return
		parameters = [float(f) * kJmol for f in fields[5:10]]
	else:
		raise RuntimeError('Unsupported dihedral formula: %s' % fields[4])
	ffdata.add('dihedral', type, 'Dihedral Coeffs', parameters)
