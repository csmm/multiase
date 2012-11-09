from lammpsbase import LAMMPSBase, SequenceType, FFData
import typing


class ClayFF(LAMMPSBase):
	
	def __init__(self, label='clayff', pair_cutoff=(10.0, 10.0), **kwargs):
		LAMMPSBase.__init__(self, label, **kwargs)
		
		try:
			ljcut, coulcut = pair_cutoff
		except:
			ljcut = pair_cutoff
			coulcut = pair_cutoff
		
		self.parameters.units          = 'real'
		self.parameters.pair_style     = 'lj/cut/coul/long %f %f' % (ljcut, coulcut)
		#self.parameters.pair_style     = 'lj/cut/coul/cut %f %f' % (ljcut, coulcut)
		self.parameters.pair_modify    = 'mix arithmetic'
		self.parameters.bond_style     = 'harmonic'
		self.parameters.angle_style    = 'harmonic'
		self.parameters.kspace_style   = 'ewald/disp 0.0001'
		
		self.ff_data = ff_data
		self.type_resolver = typing.TypeResolver(clayff_types)
		
	def atom_types(self, atoms):
		self.type_resolver.resolve_atoms(atoms)
		return atoms.info['atom_types']
	
	def set_charges(self, atoms, atom_types):
		atoms.set_charges([charges[tp] for tp in atom_types])

#sym  q      D0        R0
nonbond_data = '''
h*    0.4100 0         1      water hydrogen
ho    0.4250 0         1      hydroxyl hydrogen
o*   -0.8200 0.1554    3.5532 water oxygen 
oh   -0.9500 0.1554    3.5532 hydroxyl oxygen 
ob   -1.0500 0.1554    3.5532 bridging oxygen 
obos -1.1808 0.1554    3.5532 bridging oxygen with octahedral substitution
obts -1.1688 0.1554    3.5532 bridging oxygen with tetrahedral substitution
obss -1.2996 0.1554    3.5532 bridging oxygen with double substitutio
ohs  -1.0808 0.1554    3.5532 hydroxyl oxygen with substitution
st    2.1000 1.8405e-6 3.7064 tetrahedral silicon 
ao    1.5750 1.3298e-6 4.7943 octahedral aluminum 
at    1.5750 1.8405e-6 3.7064 tetrahedral aluminum
mgo   1.3600 9.0298e-7 5.9090 octahedral magnesium
mgh   1.0500 9.0298e-7 5.9090 hydroxide magnesium
cao   1.3600 5.0298e-6 6.2484 octahedral calcium
cah   1.0500 5.0298e-6 6.2428 hydroxide calcium
feo   1.5750 9.0298e-6 5.5070 octahedral iron
lio   0.5250 9.0298e-6 4.7257 octahedral lithium
Na    1.0 0.1301       2.6378 aqueous sodium ion
K     1.0 0.1000       3.7423 aqueous potassium ion
Cs    1.0 0.1000       4.3002 aqueous cesium ion
Ca    2.0 0.1000       3.2237 aqueous cesium ion
Ba    2.0 0.0470       4.2840 aqueous barium ion
Cl   -1.0 0.1001       4.9388 aqueous chloride ion
'''

#      k1       r0
bond_stretch_data = '''
o*  h* 554.1349 1.0000
oh  ho 554.1349 1.0000
ohs ho 554.1349 1.0000
'''

#        k2         theta0
angle_bend_data = '''
h*  o*   h* 45.7696 109.47
ao  oh   ho 30.0    109.47
at  oh   ho 30.0    109.47
mgo oh   ho 30.0    109.47
mhg oh   ho 30.0    109.47
cao oh   ho 30.0    109.47
cah oh   ho 30.0    109.47
feo oh   ho 30.0    109.47
lio oh   ho 30.0    109.47
ao  ohs  ho 30.0    109.47
at  ohs  ho 30.0    109.47
mgo ohs  ho 30.0    109.47
mhg ohs  ho 30.0    109.47
cao ohs  ho 30.0    109.47
cah ohs  ho 30.0    109.47
feo ohs  ho 30.0    109.47
lio ohs  ho 30.0    109.47
'''

ff_data = FFData()

rho_fac = 2**(-1./6)
charges = {}
for line in nonbond_data.strip().split('\n'):
	fields = line.split()
	charges[fields[0]] = float(fields[1])
	parameters = [float(fields[2]), float(fields[3])*rho_fac]
	ff_data.add('atom', fields[0], 'Pair Coeffs', parameters)

for line in bond_stretch_data.strip().split('\n'):
	fields = line.split()
	parameters = map(float, fields[2:4])
	ff_data.add('bond', SequenceType(fields[:2]), 'Bond Coeffs', parameters)

for line in angle_bend_data.strip().split('\n'):
	fields = line.split()
	parameters = map(float, fields[3:5])
	ff_data.add('angle', SequenceType(fields[:3]), 'Angle Coeffs', parameters)

clayff_types = r"""

type: h*
	! Water hydrogen
	template: [H [O[H]]]

type: ho
	! Hydroxyl hydrogen
	template: [H (O(^H))]
	
type: o*
	! Water oxygen
	template: [O [H][H]]
	
type: oh
	! Hydroxyl oxygen
	template: (O [H](^H))
	
type: ob
	! Bridging oxygen
	template: (O (^H)(^H))
	
type: obos
	! Bridging oxygen with octahedral substitution
	template: [O (Si)(Al)(^SiAl(O(H)))]
	
type: obts
	! Bridging oxygen with tetrahedral substitution
	template: (O (Si)(^SiAl))

type: obss
	! Bridging oxygen with double substitution
	template: (O (SiAl)(^SiAl)(^SiAl))
	
type: ohs
	! Hydroxyl oxygen with substitution
	template: (O [H](^Al))

# Oxygen precedences
precedence: ((o*) (oh (ohs)) (ob (obos) (obts) (obss)))
	
type: st
	! Tetrahedral silicon 
	template: [Si (*)(*)(*)(*)]
	
type: ao
	! Octahedral aluminum 
	template: [Al (*)(*)(*)(*)(*)(*)]
	
type: at
	! Tetrahedral aluminum
	template: [Al (*)(*)(*)(*)]
	
type: mgo
	! Octahedral magnesium
	template: [Mg (*)(*)(*)(*)(*)(*)]

# I don't know what this is
#type: mgh
#	! hydroxide magnesium

type: cao
	! Octahedral calcium
	template: [Ca (*)(*)(*)(*)(*)(*)]
	
#type: cah
#	! Hydroxide calcium

type: feo
	! Octahedral iron
	template: [Fe (*)(*)(*)(*)(*)(*)]
	
type: lio
	! Octahedral lithium
	template: [Li (*)(*)(*)(*)(*)(*)]
	
type: Na
	! Aqueous sodium ion
	template: (Na)
	
type: K
	! Aqueous potassium ion
	template: (K)
	
type: Cs
	! Aqueous cesium ion
	template: (Cs)

type: Ca
	! Aqueous calcium ion
	template: (Ca)

type: Ba
	! Aqueous barium ion
	template: (Ba)
	
type: Cl
	! Aqueous chloride ion
	template: (Cl)
"""
