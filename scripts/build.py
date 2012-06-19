from __future__ import division

import numpy as np
from numpy import *
random.seed(None)
import exportCML
import exportCOMPASS
from dataTypes import Atom, Bond

# Constants
cc_vector = array((.154, 0, 0)) #nm
c_angle = arccos(-1/3)

ch_length = .1094
co_length = .143
oh_length = .097
# Rotation matrix generation

def xRotation(th):
	c = cos(th)
	s = sin(th)
	return array([[1,0,0], [0,c,-s], [0,s,c]])
	
def yRotation(th):
	c = cos(th)
	s = sin(th)
	return array([[c,0,s], [0,1,0], [-s,0,c]])
	
def zRotation(th):
	c = cos(th)
	s = sin(th)
	return array([[c,-s,0], [s,c,0], [0,0,1]])


class amorphousBuilder(object):
	
	def __init__(self, boxDimensions):
		self.box = array(boxDimensions)
		self.atoms = []
		self.bonds = []
		self.temperature = 300

	def createCH3(self):
		carbon = Atom('C', zeros((3)))
		h_coords = array((-cos(c_angle), 0, sin(c_angle)))*ch_length
		h1 = Atom('H', h_coords)
		h2 = h1.copy()
		rotation = xRotation(2/3*pi)
		h2.transform(rotation)
		h3 = h2.copy()
		h3.transform(rotation)
		return [carbon, h1, h2, h3], [Bond(carbon, h1), Bond(carbon, h2), Bond(carbon, h3)]
		

	def createCH2(self):
		carbon = Atom('C', zeros((3)))
		h_coords = array((-cos(c_angle), 0, sin(c_angle)))*ch_length
		h1 = Atom('H', h_coords)
		rotation = xRotation(2/3*pi)
		h1.transform(rotation)
		h2 = h1.copy()
		h2.transform(rotation)
		return [carbon, h1, h2], [Bond(carbon, h1), Bond(carbon, h2)]

		
	def PolyEthyleneMonomerGenerator(self, N):
		for i in range(N):
			yield self.createCH2()


	def PVAMonomerGenerator(self, N):
		rotation = xRotation(2/3*pi)
		transVec = array((-cos(c_angle), 0, sin(c_angle)))
		
		# Note: one repeat unit consist of 2 carbons ("monomers")
		for i in range(N*2-1):
			if i % 2 == 1: 
				yield self.createCH2()
				continue
			
			handedness = (random.randint(2)*2-1) 
			rotation = xRotation(handedness*2/3*pi)
		
			h1_coords = dot(rotation.T, transVec*ch_length)
			
			o_coords = dot(rotation, transVec*co_length)
			h2_coords = dot(rotation, transVec*co_length+array((oh_length,0,0)))
			
			carbon = Atom('C', zeros((3)))
			h1 = Atom('H', h1_coords)
			oxygen = Atom('O', o_coords)
			h2 = Atom('H', h2_coords)
			yield [carbon, h1, oxygen, h2], [Bond(carbon, h1), Bond(carbon, oxygen), Bond(oxygen, h2)]


	# Randomly choose the conformation based on Boltzmann distribution
	def chooseConformation(self, energies):
		kT = self.temperature * 1.3806503e-23
		
		exponents = -array(energies)/kT
		
		# Avoid overflows by using logsumexp
		from scipy.maxentropy.maxentutils import logsumexp
		logsum = logsumexp(exponents)
		P = exp(exponents-logsum)
		
		print (array(energies)-energies[0])*6.0221415e23
		print P
		rval = random.uniform()
		if rval < P[0]: return 0
		if rval > 1-P[2]: return 2
		return 1
	

	def getPairInteraction(self, newAtoms):
		from lammps import lammps
		
		lmp = lammps()
		lmp.command('units real')
		lmp.command('region mybox block 0 {0} 0 {1} 0 {2} units box'.format(*self.box*10))
		lmp.command('create_box 3 mybox')
		
		types = dict(C=1, H=2, O=3)
		for atom in self.atoms + newAtoms:
			type = types[atom.type]
			lmp.command('create_atoms {0} single {1} {2} {3} units box'.format(type, *atom.pos*10))
			
		lmp.command('group system id <= {0}'.format(len(self.atoms)))
		lmp.command('group newAtoms id > {0}'.format(len(self.atoms)))
		
		lmp.command('pair_style lj/class2 10')
		lmp.command('pair_coeff 1 1 0.0620 3.8540')
		lmp.command('pair_coeff 2 2 0.0230 2.878')
		lmp.command('pair_coeff 3 3 0.0800 3.3')
		
		lmp.command('mass * 1')	
		lmp.command('neighbor 1 bin')
		lmp.command('atom_modify sort 0 0')
		lmp.command('compute energy newAtoms group/group system')
		lmp.command('thermo_style custom c_energy')
		lmp.command('run 0')
		conversionFactor = 6.9476946e-21
		energy = conversionFactor * lmp.extract_compute('energy', 0, 0)
		
		return energy
	
	
	# Main function
	def polymerGenerator(self, monomerGenerator):
		
		# Transformation from working coordinates to real coordinates
		# rotation - rotation matrix
		# translation - vector
		# Randomize initial values
		tx, ty, tz = random.uniform(0, 2*pi, size=3)
		rotation = dot(xRotation(tx), dot(yRotation(ty), zRotation(tz)))
		translation = array([random.uniform(0, dim) for dim in self.box])
		def translateOneStep():
			return (translation + dot(rotation, cc_vector)) % self.box
		
		# Create the seed
		atoms, bonds = self.createCH3()
		for atom in atoms: 
			atom.transform(zRotation(pi))
			atom.transform(rotation, translation)
			atom.periodicBox(self.box)
		self.atoms += atoms
		self.bonds += bonds
		lastCarbon = self.atoms[0]
		
		# Create the chain
		for newAtoms, newBonds in monomerGenerator:
			translation = translateOneStep()
			
			# Connect the newly created monomer to the chain
			self.bonds.append(Bond(lastCarbon, newAtoms[0]))
			lastCarbon = newAtoms[0]
			
			# Calculate energies for trans, gauche+, gauche-
			conformationEnergies = []
			for conf in (0, 1, 2):
				# Make a copy of the current monomer
				testAtoms = [atom.copy() for atom in newAtoms]
				
				# Use a CH3 as a placeholder for the next monomer for energy calculations
				# Using the correct monomer unit wouldn't make sense because we don't know
				# its conformation
				placeHolder, phBonds = self.createCH3()
				phRot = yRotation(-pi+c_angle)
				for atom in placeHolder: 
					atom.transform(phRot, translation=dot(phRot, cc_vector))
				testAtoms += placeHolder
				
				# Set the conformation (0 = trans) and calculate energy
				rot = dot(rotation, xRotation(pi * (1+conf*2/3)))
				for atom in testAtoms:
					atom.transform(rot, translation)
					atom.periodicBox(self.box)
				conformationEnergies.append(self.getPairInteraction(testAtoms))
				
			conformation = self.chooseConformation(conformationEnergies)
			rotation = dot(rotation, xRotation(pi * (1+conformation*2/3)))
			for atom in newAtoms:
				atom.transform(rotation, translation)
				atom.periodicBox(self.box)
			self.atoms += newAtoms
			self.bonds += newBonds
			rotation = dot(rotation, yRotation(-pi+c_angle))
			
			# Return the energy for records
			yield conformationEnergies[conformation]
		
		# Add the terminal group
		endAtoms, endBonds = self.createCH3()
		translation = translateOneStep()
		for atom in endAtoms: 
			atom.transform(dot(rotation, xRotation(pi)), translation)
			atom.periodicBox(self.box)
		self.bonds.append(Bond(lastCarbon, endAtoms[0]))
		self.atoms += endAtoms
		self.bonds += endBonds
		

def plotEnergies(energyLists):
	import pylab
	electronVolt = 1.60217646e-19
	for en in energyLists:
		pylab.semilogy(array(en)/electronVolt, '.')
	pylab.xlabel('repeat unit')
	pylab.ylabel('Energy (eV)')
	pylab.show()
		
boxSize = 2.0 # nm
nChains = 1

avogadro = 6.0221e23
density = 1.25 *avogadro/46/1e21
#N = 50
N = int(density*boxSize**3 / nChains)

builder = amorphousBuilder(ones((3))*boxSize)

generators = []
for i in range(nChains):
	#generators.append(builder.polymerGenerator(builder.PolyEthyleneMonomerGenerator(N)))
	generators.append(builder.polymerGenerator(builder.PVAMonomerGenerator(N)))

confEnergyLists = []
for generator in generators: confEnergyLists.append([])
	
while len(generators) > 0:
	for generator, confEnergyList in zip(generators, confEnergyLists):
		try:
			confEnergyList.append(generator.next())
		except StopIteration:
			generators.remove(generator)

#plotEnergies(confEnergyLists)
#confEnergies = [energy for energy in polymerGen]
atoms, bonds = builder.atoms, builder.bonds

exportCML.write(atoms, bonds, open('testMolecule.cml', 'w'))
exportCOMPASS.write(ones((3))*boxSize, atoms, bonds, open('testMolecule.data', 'w'))

data = dict(atoms = atoms, bonds = bonds)
import pickle
pickle.dump(data, open('testMolecule.pickle', 'w'))
