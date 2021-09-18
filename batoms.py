"""Definition of the Batoms class.

This module defines the batoms object in the blase package.

"""

from ase import Atom, Atoms
from ase.data import chemical_symbols, covalent_radii
from blase.batom import Batom
from blase.bondsetting import Bondsetting
import bpy
import bmesh
from mathutils import Vector
from copy import copy
from blase.tools import get_bondpairs, get_cell_vertices, get_bond_kind, \
                        get_polyhedra_kind, search_pbc, get_bbox
from blase.bdraw import draw_cell, draw_bond_kind, draw_polyhedra_kind, draw_text, draw_isosurface, bond_source, cylinder_mesh_from_instance, clean_default
from blase.btools import object_mode
import numpy as np
import time

subcollections = ['atom', 'bond', 'instancer', 'instancer_atom', 'cell', 'polyhedra', 'isosurface', 'virtual', 'boundary', 'text']

class Batoms():
    """Batoms Class

    In Blender, the Batoms object and collections are organised in the following way, 
    take water molecule as a example:

    * h2o                            # main collection

      * h2o_atom                   # atoms collection
    
        * atom_h2o_H                  # atoms object
    
        * atom_h2o_O                  # atoms object

      * h2o_bond                    # bonds collection

        * bond_h2o_H                 --bond object
    
        * bond_h2o_O                 --bond object
    
      * h2o_cell                     --cell collection
    
      * h2o_instancer              --instancer collection
    
        * instancer_atom_h2o_H               --sphere object
    
        * instancer_atom_h2o_O               --sphere object
    
    Then, a Batoms object is linked to this main collection in Blender. 

    Parameters:

    atoms: ase.atoms.Atoms object
        The atomic structure built using ASE
    coll: bpy.types.Collection
        The main Blender collection. 
        This collection has a ``blase`` property, which has five parameters:
            * is_batoms
            * pbc
            * unit cell
            * boundary
            * model_type
    model_type: str
        enum in '0', '1', '2', '3'
    pbc: Bool or list of Bool
        Periodic boundary conditions
    cell: array
        Unit cell size
    boundary:  list 
        search atoms at the boundary
    add_bonds: dict
        add bonds not in the default
    
    Examples:
    >>> from blase.batoms import Batoms
    >>> co = Batoms(['C', 'O'], [[[0, 0, 0]], [1.2, 0, 0]])
    >>> from ase.build import molecule
    >>> from blase.batoms import Batoms
    >>> h2o = molecule('H2O')
    >>> h2o = Batoms(h2o)
    >>> h2o.draw()
    """
    

    def __init__(self, 
                species_dict = {},
                atoms = None, 
                structure = None, 
                from_collection = None,
                label = None,
                model_type = '0', 
                boundary = [0.0, 0.0, 0.0],
                bondlist = None,
                add_bonds = {}, 
                remove_bonds = {},
                hydrogen_bond = None,
                polyhedra_dict = {},
                show_unit_cell = True,
                isosurface = [],
                kind_props = {},
                color_style = 'JMOL',
                material_style = 'blase',
                bsdf_inputs = None,
                movie = False,
                draw = True, 
                 ):
        #
        self.batoms_boundary = {}
        self.batoms_bond = {}
        self.scene = bpy.context.scene
        self.bondlist = bondlist
        self.add_bonds = add_bonds
        self.remove_bonds = remove_bonds
        self.hydrogen_bond = hydrogen_bond
        self.polyhedra_dict = polyhedra_dict
        self.isosurface = isosurface
        self.kind_props = kind_props
        self.label = label
        self.color_style = color_style
        self.material_style = material_style
        self.bsdf_inputs = bsdf_inputs
        self.bondsetting = Bondsetting(self.label)
        if species_dict:
            self.set_collection(model_type, boundary)
            self.build_batoms(species_dict)
            bondtable = self.get_bondtable(cutoff=1.1, add_bonds=add_bonds, remove_bonds=remove_bonds)
            for key, value in bondtable.items():
                self.bondsetting[key] = value
        elif atoms:
            if isinstance(atoms, list):
                self.images = atoms
                atoms = self.images[0]
            if not self.label:
                self.label = atoms.symbols.formula.format('abc')
            self.set_collection(model_type, boundary)
            self.from_ase(atoms)
            bondtable = self.get_bondtable(cutoff=1.0, add_bonds=add_bonds, remove_bonds=remove_bonds)
            for key, value in bondtable.items():
                self.bondsetting[key] = value
        elif structure:
            if isinstance(structure, list):
                self.images = structure
                structure = self.images[0]
            # if not self.label:
                # self.label = atoms.symbols.formula.format('abc')
            self.set_collection(model_type, boundary)
            self.from_pymatgen(structure)
            bondtable = self.get_bondtable(cutoff=1.0, add_bonds=add_bonds, remove_bonds=remove_bonds)
            for key, value in bondtable.items():
                self.bondsetting[key] = value
        elif from_collection:
            print('Build from collection')
            self.from_collection(from_collection)
        else:
            raise Exception("Failed, species_dict, atoms or coll should be provided!"%self.label)
        self.coll.blase.show_unit_cell = show_unit_cell
        if draw:
            self.draw()
        if movie:
            self.load_frames()
    
    
    def build_batoms(self, species_dict):
        """
        """
        for species, data in species_dict.items():
            if species not in self.kind_props: self.kind_props[species] = {}
            if isinstance(data, list):
                ba = Batom(self.label, species, data, props = self.kind_props[species], material_style=self.material_style, bsdf_inputs=self.bsdf_inputs, color_style=self.color_style)
                self.coll.children['%s_atom'%self.label].objects.link(data.batom)
                self.coll.children['%s_instancer'%self.label].objects.link(data.instancer)
            elif isinstance(data, Batom):
                self.coll.children['%s_atom'%self.label].objects.link(data.batom)
                self.coll.children['%s_instancer'%self.label].objects.link(data.instancer)
    def from_ase(self, atoms):
        """
        Import structure from ASE atoms.
        """
        if 'species' not in atoms.info:
            atoms.info['species'] = atoms.get_chemical_symbols()
        species_list = list(set(atoms.info['species']))
        for species in species_list:
            indices = [index for index, x in enumerate(atoms.info['species']) if x == species]
            if species not in self.kind_props: self.kind_props[species] = {}
            ba = Batom(self.label, species, atoms.positions[indices], props = self.kind_props[species], material_style=self.material_style, bsdf_inputs=self.bsdf_inputs, color_style=self.color_style)
            self.coll.children['%s_atom'%self.label].objects.link(ba.batom)
            self.coll.children['%s_instancer'%self.label].objects.link(ba.instancer)
        self.coll.is_batoms = True
        self.coll.blase.pbc = self.npbool2bool(atoms.pbc)
        self.coll.blase.cell = atoms.cell[:].flatten()
    def from_pymatgen(self, structure):
        """
        Import structure from Pymatgen structure.
        """
        symbols = [str(site.specie.symbol) for site in structure]
        if hasattr(structure, "lattice"):
            cell = structure.lattice.matrix
            pbc = True
        else:
            cell = None
            pbc = None
        species_list = list(set(symbols))
        for species in species_list:
            positions = [structure[index].coods for index, x in enumerate(symbols) if x == species]
            ba = Batom(self.label, species, positions, material_style=self.material_style, bsdf_inputs=self.bsdf_inputs, color_style=self.color_style)
            self.coll.children['%s_atom'%self.label].objects.link(ba.batom)
            self.coll.children['%s_instancer'%self.label].objects.link(ba.instancer)
        self.coll.is_batoms = True
        self.coll.blase.pbc = self.npbool2bool(pbc)
        self.coll.blase.cell = cell[:].flatten()
    def from_collection(self, collection_name):
        """
        """
        if collection_name not in bpy.data.collections:
            raise Exception("%s is not a collection!"%collection_name)
        elif not bpy.data.collections[collection_name].is_batoms:
            raise Exception("%s is not Batoms collection!"%collection_name)
        self.label = collection_name
    def from_pymatgen(self):
        """
        """
    def npbool2bool(self, pbc):
        """
        """
        newpbc = []
        for i in range(3):
            if pbc[i]:
                newpbc.append(True)
            else:
                newpbc.append(False)
        return newpbc
    def coll2atoms(self):
        """
        build ASE Atoms from a blase's main collection
        """
        atoms = Atoms()
        atoms.info['species'] = []
        coll = self.coll
        # atoms
        for obj in coll.children['%s_atom'%coll.name].all_objects:
            ele = obj.name.split('_')[2]
            species = '_'.join(obj.name.split('_')[2:])
            if len(obj.children) != 0:
                for vertex in obj.data.vertices:
                    location = obj.matrix_world @ vertex.co
                    atoms.append(Atom(symbol = ele, position = location))
                    atoms.info['species'].append(species)
            else:
                if not obj.parent:
                    location = obj.location
                    atoms.append(Atom(ele, location))
                    atoms.info['species'].append(species)
        # cell
        coll_cell = coll.children['%s_cell'%coll.name]
        cell_vertexs = []
        if 'cell_%s_point'%self.label in coll_cell.all_objects.keys():
            obj = coll_cell.all_objects['cell_%s_point'%self.label]
            for vertex in obj.data.vertices:
                location = obj.matrix_world @ vertex.co
                cell_vertexs.append(location)
        if cell_vertexs:
            cell = [cell_vertexs[4], cell_vertexs[2], cell_vertexs[1]]
            atoms.cell = cell
            atoms.pbc = coll.blase.pbc
    def set_collection(self, model_type = '0', boundary = [0, 0, 0]):
        """
        build main collection and its child collections.
        """
        for coll in bpy.data.collections:
            if self.label == coll.name:
                raise Exception("Failed, the name %s already in use!"%self.label)
        coll = bpy.data.collections.new(self.label)
        self.coll.blase.is_batoms = True
        self.scene.collection.children.link(self.coll)
        for sub_name in subcollections:
            subcoll = bpy.data.collections.new('%s_%s'%(self.label, sub_name))
            self.coll.children.link(subcoll)
        self.coll.blase.model_type = model_type
        self.coll.blase.boundary = boundary
    
    def draw_cell(self):
        """
        Draw unit cell
        """
        object_mode()
        cell_vertices = get_cell_vertices(self.coll.blase.cell)
        if np.max(abs(cell_vertices)) < 1e-6:
            return 0
        print('--------------Draw cell--------------')
        self.clean_blase_objects('cell')
        for obj in self.coll.children['%s_cell'%self.label].all_objects:
            bpy.data.objects.remove(obj)
        if self.show_unit_cell:
            draw_cell(self.coll.children['%s_cell'%self.label], cell_vertices, label = self.label)
    def draw_bonds(self):
        """
        Draw bonds.

        Parameters:

        cutoff: float
            cutoff used to build bond pairs.
        """
        print('--------------Draw bonds--------------')
        # if not self.bondlist:
        object_mode()
        atoms = self.get_atoms_boundary()
        self.bondlist = get_bondpairs(atoms, self.bondsetting.data)
        if self.hydrogen_bond:
            self.hydrogen_bondlist = get_bondpairs(self.atoms, cutoff = {('O', 'H'): self.hydrogen_bond})
        self.calc_bond_data(atoms, self.bondlist)
        for species, bond_data in self.bond_kinds.items():
            print('Bond %s'%species)
            draw_bond_kind(species, bond_data, label = self.label, 
                        coll = self.coll.children['%s_bond'%self.label])
    def draw_polyhedras(self):
        """
        Draw bonds.
        Parameters:

        cutoff: float
            cutoff used to build bond pairs.
        """
        object_mode()
        print('--------------Draw polyhedras--------------')
        self.calc_polyhedra_data(atoms = self.atoms, bondlist = self.bondlist)
        for species, polyhedra_data in self.polyhedra_kinds.items():
            print('Polyhedra %s'%species)
            draw_polyhedra_kind(species, polyhedra_data, label = self.label,
                        coll = self.coll.children['%s_polyhedra'%self.label])
    def draw_isosurface(self, isosurface = []):
        """
        Draw bonds.

        Parameters:

        isosurface: list
            isosurface data.
        """
        object_mode()
        if not isosurface:
            isosurface = self.isosurface
        volume = self.isosurface[0]
        icolor = 0
        if len(self.isosurface) == 1:
            draw_isosurface(self.coll.children['%s_isosurface'%self.label], volume, cell = self.atoms.cell, level=None, icolor = icolor)
        for level in self.isosurface[1:]:
            draw_isosurface(self.coll.children['%s_isosurface'%self.label], volume, cell = self.atoms.cell, level=level, icolor = icolor)
            icolor += 1
    def draw_cavity(self, radius):
        """
        cavity
        for porous materials
        >>> from ase.io import read
        >>> atoms = read('docs/source/_static/datas/mof-5.cif')
        """
        from blase.tools import find_cage
        object_mode()
        self.clean_blase_objects('virtual')
        positions = find_cage(self.cell, self.atoms.positions, radius)
        ba = Batom(self.label, 'Au_cavity', positions, scale = radius/2.8, material_style='blase', bsdf_inputs=self.bsdf_inputs, color_style=self.color_style)
        self.coll.children['%s_virtual'%self.label].objects.link(ba.batom)
        self.coll.children['%s_virtual'%self.label].objects.link(ba.instancer)
    def clean_blase_objects(self, object):
        """
        remove all bond object in the bond collection
        """
        for obj in self.coll.children['%s_%s'%(self.label, object)].all_objects:
            bpy.data.objects.remove(obj)
    def show_index(self, index_type = 0):
        """
        """
        bpy.context.preferences.view.show_developer_ui = True
        for a in bpy.context.screen.areas:
            if a.type == 'VIEW_3D':
                overlay = a.spaces.active.overlay
                overlay.show_extra_indices = True
                
    def draw(self, model_type = None):
        """
        Draw atoms, bonds, polyhedra.

        Parameters:

        model_type: str

        """
        if not model_type:
            model_type = self.model_type
        else:
            self.model_type = model_type
        self.draw_cell()
        bpy.ops.ed.undo_push()
        self.clean_blase_objects('bond')
        bpy.ops.ed.undo_push()
        self.clean_blase_objects('polyhedra')
        bpy.ops.ed.undo_push()
        if model_type == '0':
            for batom in self.batoms.values():
                batom.scale = 1.0
            for batom in self.batoms_boundary.values():
                batom.scale = 1.0
            for batom in self.batoms_bond.values():
                batom.scale = 1.0
        elif model_type == '1':
            # view(images)
            for batom in self.batoms.values():
                batom.scale = 0.4
            for batom in self.batoms_boundary.values():
                batom.scale = 0.4
            for batom in self.batoms_bond.values():
                batom.scale = 0.4
            self.draw_bonds()
        elif model_type == '2':
            for batom in self.batoms.values():
                batom.scale = 0.4
            for batom in self.batoms_boundary.values():
                batom.scale = 0.4
            for batom in self.batoms_bond.values():
                batom.scale = 0.4
            self.draw_bonds()
            self.draw_polyhedras()
        elif model_type == '3':
            for batom in self.batoms.values():
                batom.scale = 0.01
            for batom in self.batoms_boundary.values():
                batom.scale = 0.01
            for batom in self.batoms_bond.values():
                batom.scale = 0.01
            # self.clean_atoms
            self.draw_bonds()
        if self.isosurface:
            self.draw_isosurface()
    def replace(self, species1, species2, index = []):
        """
        replace atoms.

        Parameters:
        

        index: list
            index of atoms will be replaced.

        species1: str

        species2: str
            atoms will be changed to this element.

        >>> from ase.build import molecule, fcc111
        >>> from blase.batoms import Batoms
        >>> pt111 = fcc111('Pt', (5, 5, 4), vacuum = 5.0)
        >>> pt111 = Batoms(atoms = pt111, label = 'pt111')
        >>> pt111.replace('Pt', 'Au', [93])
        >>> pt111.replace('Pt', 'Au', range(20))

        """
        # if kind exists, merger, otherwise build a new kind and add.
        object_mode()
        name = 'atom_%s_%s'%(self.label, species2)
        positions = [self.batoms[species1][i] for i in index]
        if species2 in self.batoms:
            self.batoms[species2].add_vertices(positions)
        else:
            ba = Batom(self.label, species2, positions, material_style=self.material_style, bsdf_inputs=self.bsdf_inputs, color_style=self.color_style)
            self.coll.children['%s_atom'%self.label].objects.link(ba.batom)
            self.coll.children['%s_instancer'%self.label].objects.link(ba.instancer)
        self.batoms[species1].delete(index)
            
    
    def delete(self, species, index = []):
        """
        delete atoms.

        species: str

        index: list
            index of atoms to be delete
        
        For example, delete the second atom in H species.
        Please note that index start from 0.

        >>> h2o.delete([1])

        """
        self.batoms[species].delete(index)
    def translate(self, displacement):
        """Translate atomic positions.

        The displacement argument is an xyz vector.

        For example, move h2o molecule by a vector [0, 0, 5]

        >>> h2o.translate([0, 0, 5])

        """
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        for obj in self.coll.all_objects:
            if 'instancer' not in obj.name:
                obj.select_set(True)
        bpy.ops.transform.translate(value=displacement)
    def rotate(self, angle, axis = 'Z', orient_type = 'GLOBAL'):
        """Rotate atomic based on a axis and an angle.

        Parameters:

        angle: float
            Angle that the atoms is rotated around the axis.
        axis: str
            'X', 'Y' or 'Z'.

        For example, rotate h2o molecule 90 degree around 'Z' axis:
        
        >>> h2o.rotate(90, 'Z')

        """
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        for coll in self.coll.children:
            for obj in coll.objects:
                obj.select_set(True)
        bpy.ops.transform.rotate(value=angle, orient_axis=axis.upper(), orient_type = orient_type)
    
    def __getitem__(self, species):
        """Return a subset of the Batom.

        species -- str, describing which batom to return.
        """

        if isinstance(species, str):
            if species not in self.batoms:
                raise SystemExit('%s is not in this structure'%species)
            return self.batoms[species]
        elif isinstance(species, list):
            raise SystemExit('dict not supported yet!')
            
    def __add__(self, other):
        self += other
        return self
    def __iadd__(self, other):
        self.extend(other)
        return self
    def __repr__(self) -> str:
        text = []
        text.append('label={0}, '.format(self.label))
        text.append('species='.format(self.cell))
        text.append('%s '%(list(self.batoms)))
        text.append('cell={0}, '.format(self.cell))
        text.append('pbc={0}'.format(self.pbc))
        text = "".join(text)
        text = "Batom(%s)"%text
        return text
    def extend(self, other):
        """
        Extend batoms object by appending batoms from *other*.
        
        >>> from ase.build import molecule, fcc111
        >>> from blase.batoms import Batoms
        >>> import numpy as np
        >>> co = molecule('CO')
        >>> co = Batoms(atoms = co, draw = True)
        >>> au = fcc111('Au', (5, 5, 4), vacuum=5.0)
        >>> au = Batoms(atoms = au, draw = True)
        >>> co.translate(au.atoms[-1].position + np.array([0, 0, 2]))
        >>> au.extend(co)
        >>> au.atoms.save('au111-co.cif')
        
        or,

        >>> au = au + co

        """
        for species, batom in other.batoms.items():
            if species in self.species:
                self.batoms[species].extend(batom)
            else:
                ba = batom.copy(self.label, species)
                self.coll.children['%s_atom'%self.label].objects.link(ba.batom)
        self.remove_collection(other.label)

    def remove_collection(self, name):
        collection = bpy.data.collections.get(name)
        for obj in collection.all_objects:
            bpy.data.objects.remove(obj, do_unlink=True)
        for coll in collection.children:
            bpy.data.collections.remove(coll)
        bpy.data.collections.remove(collection)
    def __imul__(self, m):
        """
        """
        for species, batom in self.batoms.items():
            batom.repeat(m, self.cell)
        self.cell = np.array([m[c] * self.cell[c] for c in range(3)])
        self.set_cell(self.cell)
        return self
    def repeat(self, rep):
        """
        Create new repeated atoms object.

        >>> from ase.build import bulk
        >>> from blase.batoms import Batoms
        >>> au = bulk('Au', cubic = True)
        >>> au = Batoms(atoms = au)
        >>> au.draw()
        >>> au.repeat([2, 2, 2])
        
        """
        self.__imul__(rep)
    def __mul__(self, rep):
        self.repeat(rep)
        return self
    def copy(self, label = None):
        """
        Return a copy.

        name: str
            The name of the copy.

        For example, copy h2o molecule:
        
        >>> h2o_new = h2o.copy(label = 'h2o_new')

        """
        if not label:
            label = self.label + 'copy'
        species_dict = {x:self.batoms[x].copy(label, x) for x in self.species}
        batoms = self.__class__(species_dict = species_dict, label = label, model_type = self.coll.blase.model_type)
        batoms.translate([0, 0, 2])
        return batoms
    def write(self, filename):
        """
        Save atoms to file.

        >>> h2o.write('h2o.cif')
        
        """
        self.atoms.write(filename)
    def update(self):
        """
        """
        pass
    def update_collection(self):
        """
        """
        self.coll = bpy.data.collections[self.label]
        self.coll2atoms()
        return 0
    @property
    def coll(self):
        return self.get_coll()
    def get_coll(self):
        return bpy.data.collections[self.label]
    @property
    def coll_atom(self):
        return self.get_coll_atom()
    def get_coll_atom(self):
        return self.coll.children['%s_atom'%self.label]
    @property
    def cell(self):
        return self.get_cell()
    @cell.setter
    def cell(self, cell):
        self.set_cell(cell)
    def get_cell(self):
        """
        Get array of cell.
        """
        return np.array(self.coll.blase.cell).reshape(3, 3)
    def set_cell(self, cell, scale_atoms=False):
        """Set unit cell vectors.

        Parameters:

        Examples:

        """
        from ase.cell import Cell
        from ase.geometry.cell import complete_cell

        cell = Cell.new(cell)
        oldcell = Cell(self.get_cell())
        self.coll.blase.cell = cell[:].flatten()
        if scale_atoms:
            M = np.linalg.solve(oldcell.complete(), cell.complete())
            for ba in self.batoms.values():
                ba.set_positions(np.dot(ba.get_positions(), M))
        self.draw_cell()
    @property
    def pbc(self):
        return self.get_pbc()
    @pbc.setter
    def pbc(self, pbc):
        self.set_pbc(pbc)
    def get_pbc(self):
        return np.array(self.coll.blase.pbc)
    def set_pbc(self, pbc):
        if isinstance(pbc, bool):
            pbc = [pbc]*3
        self.coll.blase.pbc = pbc
    @property
    def boundary(self):
        return self.get_boundary()
    @boundary.setter
    def boundary(self, boundary):
        self.set_boundary(boundary)
    def get_boundary(self):
        return np.array(self.coll.blase.boundary)
    def set_boundary(self, boundary):
        """
        >>> from blase.batoms import Batoms
        >>> from ase.io import read
        >>> atoms = read('docs/source/_static/datas/tio2.cif')
        >>> tio2 = Batoms(label = 'tio2', atoms = atoms, model_type = '2', polyhedra_dict = {'Ti': ['O']}, color_style="VESTA")
        >>> tio2.boundary = 0.5
        """
        self.clean_blase_objects('boundary')
        if isinstance(boundary, (int, float)):
            boundary = [boundary]*3
        flag =  np.array(boundary[:]) > 0.0
        if self.atoms.pbc.any() and flag.any():
            for species, batom in self.batoms.items():
                positions = search_pbc(batom.positions, self.cell, boundary)
                ba = Batom(self.label, '%s_bd'%species, positions, scale = batom.scale)
                self.coll.children['%s_boundary'%self.label].objects.link(ba.batom)
                self.batoms_boundary['%s_bd'%species] = ba
        self.coll.blase.boundary = boundary
        # todo: update bond and polyhedra
    @property
    def model_type(self):
        return self.get_model_type()
    @model_type.setter
    def model_type(self, model_type):
        self.set_model_type(model_type)
    def get_model_type(self):
        return np.array(self.coll.blase.model_type)
    def set_model_type(self, model_type):
        self.coll.blase.model_type = str(model_type)
        self.draw()
    @property
    def show_unit_cell(self):
        return self.get_show_unit_cell()
    @show_unit_cell.setter
    def show_unit_cell(self, show_unit_cell):
        self.set_show_unit_cell(show_unit_cell)
    def get_show_unit_cell(self):
        return self.coll.blase.show_unit_cell
    def set_show_unit_cell(self, show_unit_cell):
        self.coll.blase.show_unit_cell = show_unit_cell
        self.draw_cell()
    @property
    def atoms(self):
        return self.get_atoms()
    def get_atoms(self):
        """
        build ASE atoms from batoms dict.
        """
        return self.batoms2atoms(self.batoms)
    def get_atoms_boundary(self):
        """
        build ASE atoms from batoms dict.
        """
        atoms = self.get_atoms()
        atoms_boundary = self.batoms2atoms(self.batoms_boundary)
        species = atoms.info['species'] + atoms_boundary.info['species']
        atoms = atoms + atoms_boundary
        atoms.info['species'] = species
        return atoms
    @property
    def species(self):
        return self.get_species()
    def get_species(self):
        """
        build species from collection.
        """
        species = []
        for ba in self.coll_atom.objects:
            species.append(ba.species)
        return species
    @property
    def batoms(self):
        return self.get_batoms()
    def get_batoms(self):
        """
        build batom dict from collection.
        """
        batoms = {}
        for ba in self.coll_atom.objects:
            batoms[ba.species] = Batom(from_batom=ba.name)
        return batoms
    def get_bondtable(self, cutoff = 1.2, add_bonds = {}, remove_bonds = {}, polyhedra_dict = {}):
        """
        """
        from blase.default_data import default_bonds
        bondtable = {}
        for species1 in self.species:
            search = False
            polyhedra = False
            radius1 = cutoff * covalent_radii[chemical_symbols.index(species1.split('_')[0])]
            for species2 in self.species:
                if species2 not in default_bonds[species1.split('_')[0]]: continue
                radius2 = cutoff * covalent_radii[chemical_symbols.index(species2.split('_')[0])]
                bondlength = radius1 + radius2
                if species1 in polyhedra_dict:
                    if species2 in polyhedra_dict[species1]:
                        polyhedra = True
                bondtable[(species1, species2)] = [bondlength, polyhedra, search]
        bondtable.update(add_bonds)
        for key in remove_bonds:
            bondtable.pop(key)
        return bondtable
    #============================================
    def get_atoms_bond(self):
        """
        build ASE atoms from batoms dict.
        """
        atoms = self.get_atoms_boundary()
        return atoms + self.batoms2atoms(self.batoms_bond)
    def batoms2atoms(self, batoms):
        object_mode()
        atoms = Atoms()
        species_list = []
        symbols = []
        positions = []
        for species, batom in batoms.items():
            if species[-3:] == '_bd': species = species[0:-3]
            if species[-3:] == '_bo': species = species[0:-3]
            species_list.extend([species]*len(batom))
            symbol = [batom.element]*len(batom)
            symbols.extend(symbol)
            positions.extend(batom.positions)
        atoms = Atoms(symbols, positions, cell = self.cell, pbc = self.pbc)
        atoms.info['species'] = species_list
        return atoms
    def draw_constraints(self):
        """
        """
        #
        constr = self.atoms.constraints
        self.constrainatoms = []
        for c in constr:
            if isinstance(c, FixAtoms):
                for n, i in enumerate(c.index):
                    self.constrainatoms += [i]
    
    def highlight_atoms(self, indexs, shape = 'sphere', radius_scale=1.4,
                           color=(0.0, 1.0, 0.0), transmit=0.6):
        """
        """
        # Draw atoms
        #
        object_mode()
        m 
        # i = 0
        for index in indexs:
            loc = self.positions[index]
            ele = self.symbols[index]
            radii = radius_scale * self.atom_kinds[ele]['radius']
            if shape == 'cube':
                bpy.ops.mesh.primitive_cube_add(location=loc, size=radii*2)
            else:
                bpy.ops.mesh.primitive_uv_sphere_add(location=loc, radius=radii)
            ball = bpy.context.view_layer.objects.active
            bpy.ops.object.shade_smooth()
            ball.data.materials.append(material)
            ball.show_transparent = True
            coll_highlight.objects.link(ball)
    def load_frames(self, images = None):
        """
        images: list
            list of atoms. All atoms show have same species and length.
        >>> from ase.io import read
        >>> from blase import Batoms
        >>> images = read('h2o-animation.xyz', index = ':')
        >>> h2o = Batoms(label = 'h2o', atoms = images)
        >>> h2o.load_frames()
        >>> h2o.render(animation = True)
        """
        if not images:
            images = self.images
        nimage = len(images)
        if len(self.atoms) != len(images[0]):
            raise Exception("Number of atoms %s is not equal to %s."%(len(self.atoms), len(images[0])))
        atoms = images[0]
        if 'species' not in atoms.info:
            atoms.info['species'] = atoms.get_chemical_symbols()
        positions = np.array([atoms.positions for atoms in images])
        for species, ba in self.batoms.items():
            index = [atom.index for atom in atoms if atoms.info['species'][atom.index] == species]
            ba.load_frames(positions[:, index])
    
    def render(self, bbox = None, output_image = None, animation = False, **kwargs):
        """
        Render the atoms, and save to a png image.

        Support all parameters for Class Blase

        >>> h2o.render(resolution_x = 1000, output_image = 'h2o.png')
        
        """
        from blase.bio import Blase
        print('--------------Render--------------')
        print('Rendering atoms')
        if not bbox:
            bbox = get_bbox(bbox = None, atoms = self.atoms)
            kwargs['bbox'] = bbox
        if not output_image:
            kwargs['output_image'] = '%s.png'%self.label
        kwargs['animation'] = animation
        # print(kwargs)
        bobj = Blase(**kwargs)
        for function in bobj.functions:
            name, paras = function
            getattr(bobj, name)(**paras)
        # bobj.load_frames()
        bobj.render()
    def calc_bond_data(self, atoms, bondlist):
        """
        """
        bond_kinds = {}
        batoms = self.batoms
        for ind1, pairs in bondlist.items():
            kind1 = atoms.info['species'][ind1]
            element = kind1.split('_')[0]
            bond_kind = get_bond_kind(element)
            if kind1 not in bond_kinds:
                bond_kinds[kind1] = bond_kind
            # print(ind1, kind, pairs)
            for bond in pairs:
                ind2, offset = bond
                kind2 = atoms.info['species'][ind2]
                R = np.dot(offset, atoms.cell)
                vec = atoms.positions[ind1] - (atoms.positions[ind2] + R)
                length = np.linalg.norm(vec)
                nvec = vec/length
                radius1 = covalent_radii[chemical_symbols.index(kind1.split('_')[0])]
                radius2 = covalent_radii[chemical_symbols.index(kind2.split('_')[0])]
                pos = [atoms.positions[ind1] - nvec*radius1*batoms[kind1].scale[0]*0.5,
                    atoms.positions[ind2] + R + nvec*radius2*batoms[kind2].scale[0]*0.5]
                center0 = (pos[0] + pos[1])/2.0
                vec = pos[0] - pos[1]
                length = np.linalg.norm(vec)
                nvec = vec/length
                nvec = nvec + 1e-8
                # verts, faces
                v1 = nvec + np.array([1.2323, 0.493749, 0.5604937284])
                v11 = v1 - np.dot(v1, nvec)*nvec
                v11 = v11/np.linalg.norm(v11)/2.828427
                v22 = np.cross(nvec, v11)*length*length
                #
                center = (center0 + pos[0])/2.0
                bond_kinds[kind1]['centers'].append(center)
                bond_kinds[kind1]['lengths'].append(length/4.0)
                bond_kinds[kind1]['normals'].append(nvec)
                nvert = len(bond_kinds[kind1]['verts'])
                bond_kinds[kind1]['verts'].append(center + v11)
                bond_kinds[kind1]['verts'].append(center - v11)
                bond_kinds[kind1]['verts'].append(center + v22)
                bond_kinds[kind1]['verts'].append(center - v22)
                bond_kinds[kind1]['faces'].append([nvert + 0, nvert + 2, nvert + 1, nvert + 3])
        for kind, bond_data in bond_kinds.items():
            self.batoms[kind].bond_data = bond_data
        self.bond_kinds = bond_kinds
    def calc_polyhedra_data(self, atoms = None, bondlist = {}, transmit = 0.8, polyhedra_dict = {}):
        """
        Two modes:
        (1) Search atoms bonded to kind
        polyhedra_dict: {'kind': ligands}
        """
        from scipy.spatial import ConvexHull
        tstart = time.time()
        polyhedra_kinds = {}
        if not polyhedra_dict:
            polyhedra_dict = {}
            for bond, data in self.bondsetting.data.items():
                if data[1]:
                    if bond[0] not in polyhedra_dict: polyhedra_dict[bond[0]] = []
                    polyhedra_dict[bond[0]].append(bond[1])
        if not atoms:
            atoms = self.atoms
        # loop center atoms
        for kind, ligand in polyhedra_dict.items():
            # print(kind, ligand)
            if kind not in polyhedra_kinds.keys():
                element = kind.split('_')[0]
                polyhedra_kinds[kind] = get_polyhedra_kind(element)
            inds = [atom.index for atom in atoms if atom.symbol == kind]
            for ind in inds:
                vertice = []
                for bond in bondlist[ind]:
                    a2, offset = bond
                    if atoms[a2].symbol in ligand:
                        temp_pos = atoms[a2].position + np.dot(offset, atoms.cell)
                        vertice.append(temp_pos)
                nverts = len(vertice)
                # print(ind, indices, nverts)
                if nverts >3:
                    # print(ind, vertice)
                    # search convex polyhedra
                    hull = ConvexHull(vertice)
                    face = hull.simplices
                    #
                    # print(ind)
                    nverts = len(polyhedra_kinds[kind]['vertices'])
                    face = face + nverts
                    edge = []
                    for f in face:
                        edge.append([f[0], f[1]])
                        edge.append([f[0], f[2]])
                        edge.append([f[1], f[2]])
                    polyhedra_kinds[kind]['vertices'] = polyhedra_kinds[kind]['vertices'] + list(vertice)
                    polyhedra_kinds[kind]['edges'] = polyhedra_kinds[kind]['edges'] + list(edge)
                    polyhedra_kinds[kind]['faces'] = polyhedra_kinds[kind]['faces'] + list(face)
                    #
                    # print('edge: ', edge)
                    for e in edge:
                        # print(e)
                        center = (polyhedra_kinds[kind]['vertices'][e[0]] + polyhedra_kinds[kind]['vertices'][e[1]])/2.0
                        vec = polyhedra_kinds[kind]['vertices'][e[0]] - polyhedra_kinds[kind]['vertices'][e[1]]
                        length = np.linalg.norm(vec)
                        nvec = vec/length
                        # print(center, nvec, length)
                        polyhedra_kinds[kind]['edge_cylinder']['lengths'].append(length/2.0)
                        polyhedra_kinds[kind]['edge_cylinder']['centers'].append(center)
                        polyhedra_kinds[kind]['edge_cylinder']['normals'].append(nvec)
        print('get_polyhedra_kind: {0:10.2f} s'.format(time.time() - tstart))
        for kind, polyhedra_data in polyhedra_kinds.items():
            self.batoms[kind].polyhedra_data = polyhedra_data
        self.polyhedra_kinds = polyhedra_kinds

    def get_distances(self, species1, i, species2, indices, mic=False):
        """
        Return distances of atom No.i with a list of atoms.

        Use mic=True to use the Minimum Image Convention.

        >>> h2o.get_distances('O', 0, 'H', [0, 1])
        """
        from ase.geometry import get_distances

        p1 = self.batoms[species1][i]
        p2 = self.batoms[species2][indices]
        cell = None
        pbc = None
        if mic:
            cell = self.cell
            pbc = self.pbc
        D, D_len = get_distances(p1, p2, cell=cell, pbc=pbc)
        D_len.shape = (-1,)
        return D_len
    def get_angle(self, species1, i1, species2, i2, species3, i3, mic=False):
        """
        Get angle in degrees between the vectors i2->i1 and
        i2->i3.
        Use mic=True to use the Minimum Image Convention and calculate the
        angle across periodic boundaries.

        >>> h2o.get_angle('H', 0, 'O', 0, 'H', 1)

        """
        from ase.geometry import get_angles
        p1 = self.batoms[species1][i1]
        p2 = self.batoms[species2][i2]
        p3 = self.batoms[species3][i3]

        v12 = p1 - p2
        v32 = p3 - p2

        cell = None
        pbc = None

        if mic:
            cell = self.cell
            pbc = self.pbc

        return get_angles([v12], [v32], cell=cell, pbc=pbc)
