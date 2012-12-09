class FFData:
	def __init__(self):
		self.atom = {}
		self.bond = {}
		self.angle = {}
		self.dihedral = {}
		self.improper = {}
		self.equivalence = {}
		self.class2 = False
		
	def add(self, group, type, title, values):
		groupdict = getattr(self, group)
		d = groupdict.setdefault(type, {})
		d[title] = values
		
	def add_equivalence(self, formal_type, **actual_types):
		self.equivalence[formal_type] = actual_types
	
	def copy(self):
		result = FFData()
		result.extend(self)
		result.class2 = self.class2
		return result
		
	def extend(self, other):
		self.atom.update(other.atom)
		self.bond.update(other.bond)
		self.angle.update(other.angle)
		self.dihedral.update(other.dihedral)
		self.improper.update(other.improper)
		self.equivalence.update(other.equivalence)

	def get(self, group, type):
		groupdict = getattr(self, group)
		try:
			return groupdict[type]
		except KeyError:
			# Test for wildcard types like * c c *
			for key in groupdict:
				if key.match(type):
					return groupdict[key]
			raise
		
	def available_tables(self, group):
		d = {}
		for entry in getattr(self, group).values():
			d.update((name, len(params)) for name, params in entry.items())
		return d.items()
	
	def get_actual_type(self, group, type):
		if not self.equivalence: return type
		try:
			if group == 'atom':
				return self.equivalence[type]['atom']
			elif type.__class__ == SequenceType:
				actuals = [self.equivalence[tp][group] for tp in type.atom_types]
				return SequenceType(actuals)
			else:
				central, others = type.get_types()
				act_central = self.equivalence[central]['improper']
				act_others = [self.equivalence[tp]['improper'] for tp in others]
				return ImproperType(central_type=act_central, other_types=act_others)
		except KeyError:
			# Silently ignore
			return type

class SequenceType:
	def __init__(self, atom_types, wildcard=None):
		if atom_types[0] < atom_types[-1]:
			self.atom_types = list(atom_types)
		else:
			self.atom_types = list(reversed(atom_types))
		if wildcard:
			for i in range(len(atom_types)):
				if self.atom_types[i] == wildcard:
					self.atom_types[i] = None
		
	def match(self, other):
		if not None in self.atom_types: return False
		match = True
		for type1, type2 in zip(self.atom_types, other.atom_types):
			if type1 != None and type1 != type2:
				match = False
		if match:
			return True
		match = True
		for type1, type2 in zip(self.atom_types, other.atom_types[::-1]):
			if type1 != None and type1 != type2:
				match = False
		return match
	def __eq__(self, other):
		return self.atom_types == other.atom_types
	def __hash__(self):
		return int(sum(hash(tp) for tp in self.atom_types) % 1000000)
	def __repr__(self):
		return repr(self.atom_types)

class ImproperType:
	def __init__(self, atom_types=None, central_type=None, other_types=None, class2=False):
		if central_type and other_types:
			self.central = central_type
			self.others = sorted(other_types)
		elif class2:
			self.central =  atom_types[1]
			self.others = sorted([atom_types[0], atom_types[2], atom_types[3]])
		else:
			self.central = atom_types[0]
			self.others = sorted(atom_types[1:])
	def get_types(self):
		return self.central, self.others
	def match(self, other):
		if not None in self.others: return False
		if self.central != other.central:
			return False
		other_copy = list(other.others)
		for tp in self.others:
			if tp != None:
				if other.others.count(tp) < self.others.count(tp):
					return False
		return True
		
	def __eq__(self, other):
		return self.central == other.central and self.others == other.others
	def __hash__(self):
		return int((hash(self.central) + sum(hash(tp) for tp in self.others)) % 100000)
	def __repr__(self):
		return '%s, %s' %(self.central, self.others)
