from __future__ import division
from scipy.constants import *

# metal style:
#mass = grams/mole
#distance = Angstroms
#time = picoseconds
#energy = eV
#velocity = Angstroms/picosecond
#force = eV/Angstrom
#torque = eV
#temperature = degrees K
#pressure = bars
#dynamic viscosity = Poise
#charge = multiple of electron charge (+1.0 is a proton)
#dipole = charge*Angstroms
#electric field = volts/Angstrom
#density = gram/cm^dim 

mole = N_A

class LAMMPSUnitConverter:
	def __init__(self, units):
		self.units = units
		
	def toASE(self, value, quantity):
		return toASE(value, quantity, self.units, 'ASE')
		
	def fromASE(self, value, quantity):
		return convert(value, quantity, 'ASE', self.units)

def convert(value, quantity, fromUnits, toUnits):
	return unitSets[fromUnits][quantity]/unitSets[toUnits][quantity] * value

def toASE(value, quantity, fromUnits):
	return convert(value, quantity, fromUnits, 'ASE')
unitSets = {}

unitSets['ASE'] = dict(
	mass = gram/mole,
	energy = eV,
	force = eV/angstrom,
	pressure = giga       #GPa
	)

unitSets['metal'] = dict(
	mass = gram/mole,
	distance = angstrom,
	time = pico,
	energy = eV,
	velocity = angstrom/pico,
	force = eV/angstrom,
	torque = eV,
	temperature = 1,
	pressure = bar,
	charge = elementary_charge,
	dipole = elementary_charge*angstrom,
	electricField = 1/angstrom,
	density = gram/centi**3
	)
		
unitSets['real'] = dict(
	mass          = gram/mole,
	distance      = angstrom,
	time          = nano,    # nanoseconds
	energy        = kilo*calorie/mole,   # kcal/mole
	velocity      = angstrom/femto,       # Angstroms/femtosecond
	force         = kilo*calorie/(mole*angstrom),
	torque        = kilo*calorie/mole,
	temperature   = 1,
	pressure      = atmosphere,
	charge        = elementary_charge,
	electricField = 1/angstrom,
	density       = gram/centi**3
	)
	

if __name__ == '__main__':
	print toMetal(15, 'energy', 'real')