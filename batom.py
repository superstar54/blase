"""Definition of the Batom class.

This module defines the Batom object in the blase package.

"""

from ase import Atom, Atoms
import bpy
import bmesh
from mathutils import Vector
from copy import copy
from blase.tools import get_atom_kind
from blase.bdraw import draw_text, draw_atom_kind, draw_bond_kind, draw_polyhedra_kind
import numpy as np
import time

subcollections = ['atom', 'bond', 'instancer', 'instancer_atom', 'polyhedra', 'isosurface', 'text']

class Vertice():
    """Vertice Class
    
    Parameters:

    species: list of str
        The atomic structure built using ASE

    positions: array

    Examples:
    >>> from blase.batom import Batom
    >>> c = Batom('C', [[0, 0, 0], [1.2, 0, 0]])
    >>> c.draw_atom()
    """
    

    def __init__(self, 
                species = None,
                position = None,
                element = None,
                x = None,
                y = None,
                z = None,
                 ):
        #
        self.species = species
        self.position = position
        self.element = element
        self.x = x
        self.y = y
        self.z = z


class Batom():
    """Batom Class
    
    Then, a Batom object is linked to this main collection in Blender. 

    Parameters:

    species: list of str
        The atomic structure built using ASE

    positions: array

    color_style: str
        "JMOL", "ASE", "VESTA"

    Examples:
    >>> from blase.batom import Batom
    >>> c = Batom('C', [[0, 0, 0], [1.2, 0, 0]])
    >>> c.draw_atom()
    """
    

    def __init__(self, 
                species,
                positions,
                name,
                scale = 1.0, 
                kind_props = {},
                color_style = 'JMOL',
                material_style = 'blase',
                bsdf_inputs = None,
                draw = False, 
                 ):
        #
        self.species = species
        self.positions = positions
        self.element = species.split('_')[0]
        self.scene = bpy.context.scene
        self.name = name
        self.kind_props = kind_props
        self.scale = scale
        self.color_style = color_style
        self.material_style = material_style
        self.bsdf_inputs = bsdf_inputs
        self.set_collection()
        self.bond_data = {}
        self.polyhedra_data = {}
        if draw:
            self.draw_atom()

    def set_collection(self):
        """
        build child collection and add it to main collections.
        """
        if self.name not in bpy.data.collections:
            self.coll = bpy.data.collections.new(self.name)
            self.coll.blase.is_blase = True
            self.scene.collection.children.link(self.coll)
            for sub_name in subcollections:
                subcoll = bpy.data.collections.new('%s_%s'%(self.name, sub_name))
                self.coll.children.link(subcoll)
        elif hasattr(bpy.data.collections[self.name], 'blase'):
            self.coll = bpy.data.collections[self.name]
        else:
            raise Exception("Failed, the name %s already in use and is not blase collection!"%self.name)
    def draw_atom(self, scale = None, kind_props = {}):
        """
        Draw atom.
        scale: float
            scale for the sphere
        kind_props: dict
            Set user defined properties for species
        """
        self.clean_blase_objects('atom')
        self.clean_blase_objects('instancer_atom')
        if scale:
            self.scale = scale
        if kind_props:
            self.kind_props = kind_props
        self.atom_kind = get_atom_kind(self.element, self.positions, scale = self.scale, props = self.kind_props, color_style = self.color_style)
        draw_atom_kind(self.species, self.coll.children['%s_atom'%self.name], 
                    self.atom_kind, name = self.name, bsdf_inputs = self.bsdf_inputs, 
                    material_style = self.material_style)
                
    def draw_bond(self, bond_data = {}):
        """
        Draw bond.
        bond_data: dict
            centers, length, normal
        """
        print('--------------Draw %s bonds--------------'%self.species)
        self.clean_blase_objects('bond')
        if bond_data:
            self.bond_data = bond_data
        draw_bond_kind(self.species, self.coll.children['%s_bond'%self.name], self.bond_data, name = self.name)
        
    def draw_polyhedra(self, polyhedra_data = {}):
        """
        Draw polyhedra.

        polyhedra_data: dict
        """
        print('--------------Draw bonds--------------')
        self.clean_blase_objects('polyhedra')
        if polyhedra_data:
            self.polyhedra_data = polyhedra_data
        if self.polyhedra_data:
            draw_polyhedra_kind(self.species, self.coll.children['%s_polyhedra'%self.name], self.polyhedra_data)
    def clean_blase_objects(self, object):
        """
        remove all bond object in the bond collection
        """
        for obj in self.coll.children['%s_%s'%(self.name, object)].all_objects:
            if obj.name == '%s_%s_%s'%(object, self.name, self.species):
                bpy.data.objects.remove(obj)
    def delete_verts(self, index = []):
        """
        delete verts
        """
        obj = self.coll.children['%s_atom'%self.name].objects['atom_%s_%s'%(self.name, self.species)]
        bm = bmesh.new()
        bm.from_mesh(obj.data)
        bm.verts.ensure_lookup_table()
        verts_select = [bm.verts[i] for i in index] 
        bmesh.ops.delete(bm, geom=verts_select, context='VERTS')
        if len(bm.verts) == 0:
            bpy.data.objects.remove(obj)
        else:
            bm.to_mesh(obj.data)
    def delete(self, index = []):
        """
        delete atom.

        index: list
            index of atoms to be delete
        
        For example, delete the second atom in h species. 
        Please note that index start from 0.

        >>> h.delete([1])

        """
        self.delete_verts(index)
    def update(self):
        """
        """
        self.draw_atom()
    
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
    
    def load_frames(self, nimages = None):
        """
        """
        # render settings
        if not nimages:
            nimages = len(self.images)
        for i in range(0, nimages):
            atom_kinds = get_atom_kinds(self.images[i])
            # bond_kinds = get_bond_kinds(self.images[i], bond_list, self.nbins)
            for kind, datas in atom_kinds.items():
                obj_atom = bpy.data.objects['atom_{0}_{1}'.format(self.name, kind)]
                nverts = len(obj_atom.data.vertices)
                for j in range(nverts):
                    obj_atom.data.vertices[j].co = datas['positions'][j]
                    obj_atom.data.vertices[j].keyframe_insert('co', frame=i + 1)
        self.scene.frame_start = 1
        self.scene.frame_end = nimages
    
    def __getitem__(self, i):
        """Return a subset of the Batom.

        i -- int, describing which atom to return.
        """

        if isinstance(i, int):
            natoms = len(self.positions)
            if i < -natoms or i >= natoms:
                raise IndexError('Index out of range.')
            vertice = {'species': self.species,
                       'element': self.element,
                       'position': self.positions[i],
                       'x': self.positions[i][0],
                       'y': self.positions[i][1],
                       'z': self.positions[i][2],
                       }
            return vertice
    def copy(self, name = None):
        """
        Return a copy.

        name: str
            The name of the copy.

        For example, copy H species:
        
        >>> h_new = h.copy(name = 'h_new')

        """
        if not name:
            name = self.name + 'copy'
        batom = self.__class__(self.species, self.positions, name)
        return batom
    def extend(self, other):
        """
        Extend batom object by appending batom from *other*.
        
        >>> from ase.build import molecule, fcc111
        >>> from blase.batoms import Batoms
        >>> import numpy as np
        >>> co = molecule('CO')
        >>> co = Batoms(atoms = co, draw = True)
        >>> co.translate([0, 0, 2])
        >>> h2o = molecule('H2O')
        >>> h2o = Batoms(atoms = h2o, draw = True)
        >>> h2o['O'].extend(co['O'])
        """
        self.positions = np.append(self.positions, other.positions, axis=0)
        self.update()
        bpy.data.objects.remove(other.coll.all_objects['atom_%s_%s'%(other.name, other.species)])
    def translate(self, displacement):
        """Translate atomic positions.

        The displacement argument is an xyz vector.

        For example, move H species molecule by a vector [0, 0, 5]

        >>> h2o["H"].translate([0, 0, 5])

        Todo only redraw atoms, not bonds
        """
        self.positions += displacement
        bpy.ops.object.select_all(action='DESELECT')
        for obj in self.coll.all_objects:
            if 'instancer' not in obj.name:
                obj.select_set(True)
        bpy.ops.transform.translate(value=displacement)