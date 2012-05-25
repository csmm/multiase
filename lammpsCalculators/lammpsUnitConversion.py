from __future__ import division
from scipy.constants import *
import ase.units

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

def convert(value, quantity, fromUnits, toUnits):
	return unitSets[fromUnits][quantity]/unitSets[toUnits][quantity] * value
	
unitSets = {}

unitSets['ASE'] = dict(
	mass = gram/mole,
	distance = angstrom,
	time = 1/ase.units.second,
	energy = eV,
	velocity = angstrom/(1/ase.units.second),
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
	time          = femto,    # femtoseconds
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
	print convert(ase.units.fs, 'time', 'ASE', 'real')