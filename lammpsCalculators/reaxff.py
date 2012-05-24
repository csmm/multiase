from lammpsBase import LAMMPSBase

class ReaxFF(LAMMPSBase):
	
	def __init__(self, label='reaxff', specorder=None, ff_file_path=None, **kwargs):
		
		LAMMPSBase.__init__(self, label, files = ['ffield.reax'], **kwargs)
		
		self.specorder = specorder
		self.parameters.atom_style = 'charge'
		self.parameters.pair_style = 'reax'
		self.parameters.units      = 'real'
		
	def prepare_calculation(self):
		self.data.clear()
		
		elements = self.atoms.get_chemical_symbols()
		species = filter(lambda s: s in elements, self.specorder)
		
		self.data.atom_types = [species.index(element)+1 for element in elements]
		
		massDict = dict((el, mass) for el, mass in zip(elements, self.atoms.get_masses()))
		self.data.masses = [massDict[el] for el in species]
		
		element_indices = [str(self.specorder.index(s)+1) for s in species]
		self.parameters.pair_coeffs = ['* * ffield.reax '+' '.join(element_indices)]
		