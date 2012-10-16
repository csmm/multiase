from __future__ import division

import numpy as np
from numpy import *
from scipy import constants
random.seed(None)
from ase.atoms import Atoms
import multiprocessing
from multiasecalc.lammps.bonds import Bonds

# Constants
cc_vector = array((1.54, 0, 0))
c_angle = arccos(-1/3)

ch_length = 1.094
co_length = 1.43
oh_length = .97
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

def addAtoms(atoms1, atoms2):
	atoms1.info['bonds'].atoms_extending(atoms2)
	atoms1.extend(atoms2)

def rotateAtoms(atoms, rotation):
	atoms.positions = dot(atoms.positions, rotation.T)
	
def pairEnergy(atoms1, atoms2, box):
	coords1 = atoms1.positions
	coords2 = atoms2.positions
	numbers1 = atoms1.numbers
	numbers2 = atoms2.numbers
	
	distances_sq = zeros((len(atoms1), len(atoms2)))
	for dim in range(3):
		X2, X1 = meshgrid(coords2[:,dim], coords1[:,dim])
		fracdist = (X2 - X1) / box[dim]
		fracdist = (fracdist+.5) % 1  - .5    # MIC
		distances_sq += (fracdist*box[dim])**2
	
	
	lj_params = {
		6: (.06, 4.0),  # C
		1:  (.04, 2.0), # H
		8: (.10, 3.4)   # O
	}
	
	sig = zeros_like(distances_sq)
	eps = ones_like(distances_sq)
	for number, values in lj_params.items():
		ind1 = where(numbers1 == number)
		ind2 = where(numbers2 == number)
		eps[ind1,:] *= values[0]
		eps[:,ind2] *= values[0]
		sig[ind1,:] += values[1]
		sig[:,ind2] += values[1]
	
	sig /= 2
	eps = np.sqrt(eps)
	
	sig_per_r = sig / np.sqrt(distances_sq)
	energies = 4*eps*(sig_per_r**12 - sig_per_r**6)
	return sum(energies)


class amorphousBuilder(object):
	
	def __init__(self, boxDimensions):
		self.box = array(boxDimensions)
		self.atoms = None
		self.temperature = 1000

	def createCH3(self):
		rot = xRotation(2/3*pi)
		mpow = np.linalg.matrix_power
		h_coords = array((-cos(c_angle), 0, sin(c_angle)))*ch_length
		pos = [np.zeros((3)), h_coords, dot(rot, h_coords), dot(mpow(rot,2), h_coords)]
		atoms = Atoms('CH3', pos)
		atoms.info['bonds'] = Bonds(atoms, pairs=((0,1), (0,2), (0,3)))
		return atoms
		

	def createCH2(self):
		rot = xRotation(2/3*pi)
		mpow = np.linalg.matrix_power
		h_coords = array((-cos(c_angle), 0, sin(c_angle)))*ch_length
		pos = [zeros((3)), dot(rot, h_coords), dot(mpow(rot,2), h_coords)]
		atoms = Atoms('CH2', pos)
		atoms.info['bonds'] = Bonds(atoms, pairs=((0,1), (0,2)))
		return atoms

		
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
			
			atoms = Atoms('CHOH', [zeros((3)), h1_coords, o_coords, h2_coords])
			atoms.info['bonds'] = Bonds(atoms, pairs=((0,1), (0,2), (2,3)))
			yield atoms


	# Randomly choose the conformation based on Boltzmann distribution
	def chooseConformation(self, energies):
		kT = self.temperature * constants.Boltzmann*constants.Avogadro/1000
		
		exponents = -array(energies)/kT
		
		# Avoid overflows by using logsumexp
		from scipy.maxentropy.maxentutils import logsumexp
		logsum = logsumexp(exponents)
		P = exp(exponents-logsum)
		
		rval = random.uniform()
		if rval < P[0]: return 0
		if rval > 1-P[2]: return 2
		return 1
	
	
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
			return (translation + dot(rotation, cc_vector))
		
		# Create the seed
		atoms = self.createCH3()
		rotateAtoms(atoms, zRotation(pi))
		rotateAtoms(atoms, rotation)
		atoms.translate(translation)
		if self.atoms:
			lastCarbon = len(self.atoms)
			addAtoms(self.atoms, atoms)
		else:
			atoms.cell = self.box
			atoms.pbc = True
			self.atoms = atoms
			lastCarbon = 0
		
		processes = multiprocessing.Pool(processes=3)
		
		# Create the chain
		for newAtoms in monomerGenerator:
			translation = translateOneStep()
		
			# Calculate energies for trans, gauche+, gauche-
			conformationEnergies = []
			energyResults = []
			for conf in (0, 1, 2):
				# Make a copy of the current monomer
				testAtoms = newAtoms.copy()
				
				# Use a CH3 as a placeholder for the next monomer for energy calculations
				# Using the correct monomer unit wouldn't make sense because we don't know
				# its conformation
				placeHolder = self.createCH3()
				phRot = yRotation(-pi+c_angle)
				rotateAtoms(placeHolder, phRot)
				placeHolder.translate(dot(phRot, cc_vector))
				testAtoms += placeHolder
				
				# Set the conformation (0 = trans) and calculate energy
				rot = dot(rotation, xRotation(pi * (1+conf*2/3)))
				rotateAtoms(testAtoms, rot)
				testAtoms.translate(translation)
				#conformationEnergies.append(pairEnergy(atoms, testAtoms, self.box))
				energyResults.append(processes.apply_async(pairEnergy, (self.atoms, testAtoms, self.box)))
			
			conformationEnergies = [result.get() for result in energyResults]
			conformation = self.chooseConformation(conformationEnergies)
			rotation = dot(rotation, xRotation(pi * (1+conformation*2/3)))
			rotateAtoms(newAtoms, rotation)
			newAtoms.translate(translation)
			
			newCarbon = len(self.atoms)
			addAtoms(self.atoms, newAtoms)
			self.atoms.info['bonds'].add(lastCarbon, newCarbon)
			lastCarbon = newCarbon
			
			rotation = dot(rotation, yRotation(-pi+c_angle))
			
			# Return the energy for records
			yield conformationEnergies[conformation]
			
		processes.close()
		processes.join()
		
		# Add the terminal group
		endAtoms = self.createCH3()
		translation = translateOneStep()
		rotateAtoms(endAtoms, dot(rotation,xRotation(pi)))
		endAtoms.translate(translation)
		
		newCarbon = len(self.atoms)
		addAtoms(self.atoms, endAtoms)
		self.atoms.info['bonds'].add(lastCarbon, newCarbon)
		
		

def plotEnergies(energyLists):
	import pylab
	electronVolt = 1.60217646e-19
	for en in energyLists:
		pylab.semilogy(array(en)/electronVolt, '.')
	pylab.xlabel('repeat unit')
	pylab.ylabel('Energy (eV)')
	pylab.show()

	
def createPVA(boxSize, nchains=1):
	box = ones((3))*boxSize
	
	angstrom = 10
	
	avogadro = 6.0221e23
	#density = 1.25 *avogadro/46/1e24
	density = 1.2 *avogadro/46/1e24
	N = int(density*boxSize**3 / nchains)
	#N = 1
	print N, 'monomers per chain'

	builder = amorphousBuilder(box)

	generators = []
	for i in range(nchains):
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

	atoms = builder.atoms
	#atoms.set_scaled_positions(atoms.get_scaled_positions())
	return atoms
	
if __name__ == '__main__':
	atoms = createPVA(12)
	from multiasecalc.lammps import COMPASS, CHARMM, dynamics
	from multiasecalc.utils import get_datafile
	from atomsview import atomsview
	from multiasecalc.lammps import compasstypes, charmmtypes
	atomsview.view(atoms, charmmtypes.data)
	atoms.calc = COMPASS(ff_file_path = get_datafile('compass.frc'), debug=True)
	#atoms.calc = CHARMM(get_datafile('par_all36_cgenff.prm'), debug=True)
	min = dynamics.LAMMPSOptimizer(atoms)
	min.run()
	from ase.visualize import view
	view(atoms)