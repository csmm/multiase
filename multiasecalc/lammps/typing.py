import re
import numpy as np

#itertools.permutations not in python 2.4

#from itertools import permutations

def permutations(iterable, r=None):
    # permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
    # permutations(range(3)) --> 012 021 102 120 201 210
    pool = tuple(iterable)
    n = len(pool)
    if r is None: r = n
    if r > n:
        return
    indices = range(n)
    cycles = range(n, n-r, -1)
    yield tuple(pool[i] for i in indices[:r])
    while n:
        for i in reversed(range(r)):
            cycles[i] -= 1
            if cycles[i] == 0:
                indices[i:] = indices[i+1:] + indices[i:i+1]
                cycles[i] = n - i
            else:
                j = cycles[i]
                indices[i], indices[-j] = indices[-j], indices[i]
                yield tuple(pool[i] for i in indices[:r])
                break
        else:
            return
            

class SyntaxError(Exception): pass
class PrecedenceError(Exception): pass
class TypingError(Exception): pass

class TypeResolver:
	def __init__(self, typedata):
		comment_exp = re.compile(r'\n\s*![^\n]\n')
		type_exp = re.compile(r'type:\s*(\S+)\s*(![\w\W]*?)?\n\s*template:\s*(.+)\n')
		precedence_exp = re.compile(r'precedence:\s*([\w\W]+?)\n\s*\n')
		
		templates = []
		for match in type_exp.finditer(typedata):
			type, doc, exp = match.groups()
			templates.append(Template(type, exp, doc))
		
		precedences = []
		for match in precedence_exp.finditer(typedata):
			precedences.append(PrecedenceTree(match.group(1)))
			
		self.templates = templates
		self.precedences = precedences
	
	def resolve(self, atom):
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
		bonded = [atoms[idx] for idx in atoms.info['bonds'][atom.index]]
		
		if parent_atom:
			bonded = [a for a in bonded if a.index != parent_atom.index]
		
		if self.closed and len(bonded) != len(self.children):
			return False
		
		if self.ring and not self.test_ring(atom):
			return False
			
		for sequence in permutations(bonded, len(self.children)):
			if np.all([ch.match(neigh, atom) for neigh, ch in zip(sequence, self.children)]):
				return True
		return False
		
	def test_ring(self, atom):
		atoms = atom.atoms
		bonds = atoms.info['bonds']
		conditions = self.ring_conditions
		
		def recursive_finder(chain, next_atom):
			
			if chain and next_atom == chain[0] and len(chain) in self.ringsize_range:
				# Test conditions for atoms except the starting atom
				for sequence in permutations(chain[1:], len(conditions)):
					if np.all([c.match(atoms[ind]) for ind, c in zip(sequence, conditions)]):
						return True
				return False
			
			if next_atom in chain:
				return False
			if len(chain) + 1 > self.ringsize_range[-1]:
				return False
				
			chain = chain + [next_atom]
			return np.any([recursive_finder(chain, atom) for atom in bonds[next_atom]])
		
		return recursive_finder([], atom.index)

		
class PrecedenceTree(Tree):
	def contains(self, template):
		return template.type == self.root or \
				np.any([ch.contains(template) for ch in self.children])
	
	def contains_all(self, templates):
		return np.all([self.contains(t) for t in templates])
	
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