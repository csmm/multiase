from PyQt4.QtCore import pyqtSlot, pyqtSignal, pyqtProperty
from PyQt4.QtCore import QAbstractListModel, QModelIndex
from PyQt4.QtCore import Qt, QObject, QUrl, QStringList
from PyQt4.QtGui import QApplication, QWidget, QColor, qRgb
from PyQt4.QtDeclarative import QDeclarativeView
import os
import numpy as np
from ase.data import covalent_radii, vdw_radii, colors, atomic_numbers
from multiasecalc.lammps.typing import TypeResolver
from multiasecalc.lammps.bonds import Bonds
from multiasecalc.lammps import charmmtypes, compasstypes
from multiasecalc import utils

qml_sourcepath = os.path.join(utils.csmm_config['CSMM_INSTALL_DIR'], 'atomsview/atomsview.qml')

def view(atoms, typedata=None):
	if typedata:
		resolver = TypeResolver(typedata)
		resolver.resolve_atoms(atoms)
	app = QApplication([])
	view = AtomsView(atoms)
	view.show()
	app.exec_()

class AtomsView(QDeclarativeView):
	def __init__(self, atoms):
		QWidget.__init__(self)
		self.view_state = ViewState(atoms)
		self.bonds_model = BondsModel(self.view_state)
		self.atoms_model = AtomsModel(self.view_state)
		
		self.rootContext().setContextProperty('viewState', self.view_state)
		self.rootContext().setContextProperty('atomsModel', self.atoms_model)
		self.rootContext().setContextProperty('bondsModel', self.bonds_model)
		self.setSource(QUrl(qml_sourcepath))
		self.setResizeMode(QDeclarativeView.SizeRootObjectToView)
		
		self.setWindowTitle('AtomsView')
		self.resize(800, 800)


class ViewState(QObject):
	updated = pyqtSignal()
	
	def __init__(self, atoms):
		QObject.__init__(self)
		self.atoms = atoms
		self.atom_coordinates = np.array(atoms.positions)
		self.translation = np.zeros(3)
		self.rotation = np.diag((1,1,1))
		self.centerAtoms()
	
	def update_coordinates(self):
		coords = np.dot(self.atoms.positions+self.translation, self.rotation.T)
		self.atom_coordinates = coords
		self.updated.emit()
	
	def centerAtoms(self):
		self.translation = -np.mean(self.atoms.positions, axis=0)
		self.update_coordinates()
		
	@pyqtSlot(float, float)
	def rotate(self, x_radians, y_radians):
		x_rot = x_rotation(y_radians)
		y_rot = y_rotation(x_radians)
		self.rotation = np.dot(np.dot(x_rot,y_rot), self.rotation)
		self.update_coordinates()
	
	@pyqtSlot(str, result='QVariantList')
	def jmolColor(self, element):
		number = atomic_numbers[str(element)]
		color = colors.jmol_colors[number]
		return [float(c) for c in color]

def x_rotation(th):
	c, s = np.cos(th), np.sin(th)
	return np.array([[1,0,0], [0,c,-s], [0,s,c]])
	
def y_rotation(th):
	c, s = np.cos(th), np.sin(th)
	return np.array([[c,0,s], [0,1,0], [-s,0,c]])


class AtomsModel(QAbstractListModel):
	
	xrole, yrole, zrole, elementrole, typerole, descriptionrole, covradiusrole = range(Qt.UserRole+1, Qt.UserRole+1+7)
	
	def __init__(self, view_state):
		QAbstractListModel.__init__(self)
		self.view_state = view_state
		view_state.updated.connect(self.changed)
		role_names = {
			self.xrole: 'atomx',
			self.yrole: 'atomy',
			self.zrole: 'atomz',
			self.elementrole: 'element',
			self.covradiusrole: 'covalentRadius',
			self.typerole: 'type',
			self.descriptionrole: 'description'
			}
		self.setRoleNames(role_names)
		
	def rowCount(self, index = QModelIndex()):
		return len(self.view_state.atoms)
		
	def data(self, index, role = Qt.DisplayRole):
		atoms = self.view_state.atoms
		row = index.row()
		if role in (self.xrole, self.yrole, self.zrole):
			return float(self.view_state.atom_coordinates[row, role - self.xrole])
		else:
			has_types = 'atom_types' in self.view_state.atoms.info
			has_doc = 'descriptions' in self.view_state.atoms.info
			return {
				self.elementrole: atoms.get_chemical_symbols()[row],
				self.typerole: atoms.info['atom_types'][row] if has_types else False,
				self.descriptionrole: atoms.info['descriptions'][row] if has_doc else False,
				self.covradiusrole: float(covalent_radii[atoms.numbers[row]])
			}[role]
		
	@pyqtSlot()
	def changed(self):
		self.dataChanged.emit(self.index(0), self.index(self.rowCount()-1))

		
class BondsModel(QAbstractListModel):
	
	x1role, y1role, z1role, x2role, y2role, z2role, element1role, element2role = range(Qt.UserRole+1, Qt.UserRole+1+8)
	
	def __init__(self, view_state):
		QAbstractListModel.__init__(self)
		self.view_state = view_state
		self.pairs = list(view_state.atoms.info['bonds'])
		view_state.updated.connect(self.changed)
		role_names = {
			self.x1role: 'x1',
			self.y1role: 'y1',
			self.z1role: 'z1',
			self.x2role: 'x2',
			self.y2role: 'y2',
			self.z2role: 'z2',
			self.element1role: 'element1',
			self.element2role: 'element2',
			}
		self.setRoleNames(role_names)
	
	def rowCount(self, index = QModelIndex()):
		return len(self.view_state.atoms.info['bonds'])
	
	def data(self, index, role = Qt.DisplayRole):
		atoms = self.view_state.atoms
		row = index.row()
		i, j = self.pairs[row]
		if role in (self.x1role, self.y1role, self.z1role):
			return float(self.view_state.atom_coordinates[i, role - self.x1role])
		if role in (self.x2role, self.y2role, self.z2role):
			return float(self.view_state.atom_coordinates[j, role - self.x2role])
		else:
			return {
				self.element1role: atoms.get_chemical_symbols()[i],
				self.element2role: atoms.get_chemical_symbols()[j],
			}[role]
	
	@pyqtSlot()
	def changed(self):
		self.dataChanged.emit(self.index(0), self.index(self.rowCount()-1))
		
if __name__ == '__main__':
	import sys
	from ase.data import s22
	
	try: 
		atoms = s22.create_s22_system(sys.argv[1])
	except:
		atoms = s22.create_s22_system('Adenine-thymine_complex_stack')
	
	try:
		if sys.argv[2] == 'charmm':
			typedata = charmmtypes.data
		if sys.argv[2] == 'compass':
			typedata = compasstypes.data
	except:
		typedata = charmmtypes.data
	
	atoms.info['bonds'] = Bonds(atoms, autodetect=True)
	type_resolver = TypeResolver(typedata)
	matches = [type_resolver.resolve(atom) for atom in atoms]
	atoms.info['atom_types'] = [m.type for m in matches]
	atoms.info['descriptions'] = [m.docstring for m in matches]
	view(atoms)