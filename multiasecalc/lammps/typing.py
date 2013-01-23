import re
import numpy as np
from bonds import Bonds

try:
	from itertools import product
except:
	#itertools.product is not in python 2.4
	def product(*args, **kwds):
		# product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
		# product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
		pools = map(tuple, args) * kwds.get('repeat', 1)
		result = [[]]
		for pool in pools:
			result = [x+[y] for x in result for y in pool]
		for prod in result:
			yield tuple(prod)
            

class SyntaxError(Exception): pass
class PrecedenceError(Exception): pass
class TypingError(Exception): pass

class TypeResolver:
	def __init__(self, typedata):
		comment_exp = re.compile(r'\n\s*![^\n]\n')
		type_exp = re.compile(r'type:\s*(\S+)\s*(![\w\W]*?)?\n\s*template:\s*(.+)\n')
		precedence_exp = re.compile(r'precedence:\s*([\w\W]+?)\n\s*\n')
		
		typedata = re.sub(r'#.*\n', '\n', typedata) # remove comment lines
		
		templates = []
		for match in type_exp.finditer(typedata):
			type, doc, exp = match.groups()
			templates.append(Template(type, exp, doc))
		
		precedences = []
		for match in precedence_exp.finditer(typedata):
			precedences.append(PrecedenceTree(match.group(1)))
			
		self.templates = templates
		self.precedences = precedences
		
		self.previous_atoms = None
	
	def resolve(self, atom):
		atoms = atom.atoms
		if atoms != self.previous_atoms or 'rings' not in atoms.info:
			atoms.info['rings'] = self.find_rings(atoms)
			self.previous_atoms = atoms.copy()
		
		matches = [temp for temp in self.templates if temp.match(atom)]
		if not matches:
			raise TypingError('Could not resolve type for atom %s' % atom)
		if len(matches) == 1: return matches[0]
		
		precs = self.precedences
		try:
			# Always use generators to make it a one-liner if possible ;) (I know, it's stupid)
			return (p.resolve(matches) for p in precs if p.contains_all(matches)).next()
		except StopIteration:
			raise PrecedenceError('%s:\nCould not resolve precedence among types %s' % (atom, matches))
			
	def resolve_atoms(self, atoms):
		templates = [self.resolve(atom) for atom in atoms]
		atoms.info['atom_types'] = [t.type for t in templates]
		atoms.info['descriptions'] = [t.docstring for t in templates]
		
	
	def find_rings(self, atoms):
		ring_size_limit = 7
		rings = []
		ring_atoms = set()
		ring_bonds = self.get_ring_bonds(atoms)
		
		def add_ring(chain):
			rings.append(chain)
			ring_atoms.update(chain)
		
		def recursive_finder(chain, next_atom):
			if len(chain) > 2 and next_atom == chain[0]:
				add_ring(chain)
				return
			
			if next_atom in chain: return
			if len(chain) + 1 > ring_size_limit: return
			
			neighbors = ring_bonds[next_atom]
			if len(neighbors) > 3: return                   # ignore rings with sp3 carbons
			
			for atom in neighbors:
				recursive_finder(chain + [next_atom], atom)
		
		for atom in range(len(atoms)):
			if not atom in ring_atoms:
				recursive_finder([], atom)
		
		return rings
		
		
	def get_ring_bonds(self, atoms):
		''' Only take bonds that are a part of a ring '''
		bm = atoms.info['bonds'].get_bond_matrix()
		prev_ends = None
		while True:
			# Remove chain end bonds until only rings remain
			sums = bm.sum(axis=0)
			ends = np.where(sums < 2)[0]
			if np.array_equal(ends, prev_ends): break
			prev_ends  = ends
			bm[:,ends] = 0
			bm[ends,:] = 0
			
		return Bonds(atoms, zip(*np.where(np.triu(bm) > 0)))


class Template:
	def __init__(self, type, template, docstring=None):
		self.type = type
		self.docstring = docstring
		self.tree = TemplateTree(template)
		
	def match(self, atom):
		return self.tree.match(atom)

	def __repr__(self): return '%s %s' % (self.type, self.docstring)

class Tree:
	def __init__(self, exp):
		exp = exp.strip()
		brackets = (exp[0], exp[-1])
		if brackets not in (('(',')'), ('[',']'), ('{','}')):
			raise SyntaxError("Invalid expression %s" % exp)
		
		positions = []
		level = 0
		for pos in range(1, len(exp)):
			if exp[pos] in ('(', '[', '{'):
				if level == 0 and exp[pos] != '{':
					positions.append(pos)
				level += 1
			elif exp[pos] in (')', ']', '}'):
				level -= 1
		
		positions.append(-1)
		self.root = exp[1:positions[0]].strip()
		children_exps = [exp[start:end] for start, end in zip(positions[:-1], positions[1:])]
		self.add_children(children_exps)
		
		self.expression = exp
		
	def add_children(self, expressions):
		self.children = [self.__class__(exp) for exp in expressions]
		
	def __str__(self):
		return self.expression


class TemplateTree(Tree):
	
	# A regular expression for an atom expression such as >CNR6{(C[-O])}
	rexp = re.compile(r'[-~=>]?(\*)?(\^)?([A-QS-Za-z]*)\s*(R)?(\d)?(\{.+\})?')
	
	# By default rings with 5-6 members are accepted
	default_ringsize_range = range(5,7)
	
	def __init__(self, exp):
		Tree.__init__(self, exp)
		self.closed = {'[': True, '(': False, '{': None}[exp[0]]
			
		try:
			m = self.rexp.match(self.root)
		except:
			raise SyntaxError('Invalid syntax: %s') % exp
		
		asterisk, exclude, elements, R, ringsize, ring_condition = m.groups()
		
		self.invert = bool(asterisk or exclude)
			
		self.elements = []
		if not asterisk:
			for m in re.finditer('[A-Z][a-z]?', elements):
				self.elements.append(m.group())
		
		self.ring = bool(R)
		if ringsize: 
			self.ringsize_range = [int(ringsize)]
		else:
			self.ringsize_range = self.default_ringsize_range
		
		if ring_condition:
			self.ring_conditions = TemplateTree(ring_condition).children
		else:
			self.ring_conditions = []
	
	def match(self, atom, parent_atom=None):
		if (atom.symbol in self.elements) == self.invert: return False
		
		atoms = atom.atoms
		bonded = atoms.info['bonds'][atom.index]
		if parent_atom: bonded.remove(parent_atom.index)
		nchildren = len(self.children)
		
		if len(bonded) < nchildren or (self.closed and len(bonded) > nchildren):
			return False
		
		if self.ring and self.test_ring(atom) == False:
			return False
			
		if nchildren == 0: return True
		
		return self.find_matching_sequence(atoms, bonded, self.children, parent_atom=atom)
		
	def test_ring(self, atom):
		atoms = atom.atoms
		conditions = self.ring_conditions
		for ring in atoms.info['rings']:
			if atom.index in ring and len(ring) in self.ringsize_range:
				if len(conditions) == 0: return True
				ring_copy = list(ring)
				ring_copy.remove(atom.index)
				# Test conditions for atoms except the current atom
				return self.find_matching_sequence(atoms, ring_copy, conditions)
		
		return False
		
	def find_matching_sequence(self, atoms, atom_indices, testers, **tester_kwargs):
		# Find a sequence in given atoms that matches the sequence of testers
		matchings = []
		for tester in testers:
			matched = [i for i in atom_indices if tester.match(atoms[i], **tester_kwargs)]
			if not matched: return False
			matchings.append(matched)
		for sequence in product(*matchings):
			if len(set(sequence)) == len(testers): return True
		return False
		
		
class PrecedenceTree(Tree):
	def contains(self, template):
		return template.type == self.root or \
				True in (ch.contains(template) for ch in self.children)
	
	def contains_all(self, templates):
		return not False in (self.contains(t) for t in templates)
	
	def resolve(self, templates):
		for ch in self.children:
			result = ch.resolve(templates)
			if result: return result
		try: 
			return (t for t in templates if t.type == self.root).next()
		except StopIteration:
			return False


if __name__ == '__main__':
	f = open('charmmtypes.template')
	resolver = TypeResolver(f)
	
	from ase.data import s22
	from bonds import Bonds

	for name in s22.s22:
		dimer = s22.create_s22_system(name)
		
		print
		print name
		dimer.info['bonds'] = Bonds(dimer, autodetect=True)
	
		for atom in dimer:
			print atom.symbol, resolver.resolve(atom)