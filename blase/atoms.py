
"""Definition of the Atoms class.

This module defines the atoms object in the blase package.


add more properties to the atoms object, including:
1) atoms kinds
2) coll
"""

from ase import Atom, Atoms
import os
import bpy
import mathutils as mu
import bmesh
from mathutils import Vector
from math import sqrt
from copy import copy
from blase.bio import Blase
from blase.btools import draw_cell, draw_atoms, draw_bonds, draw_polyhedras, draw_isosurface, bond_source, cylinder_mesh_from_instance, clean_default
import numpy as np


class Batoms(Atoms):
    """Atoms object.

    Similar to ASE atoms object

    Parameters:

    
    """

    blase_objtype = 'atoms'

    def __init__(self, symbols=None,
                 positions=None, name = 'BAtoms', **kwargs):
        Atoms.__init__(self, symbols=symbols,
                 positions=positions, **kwargs)
        self.name = name
        print(self)
        if not self.name:
            self.name = self.atoms.symbols.formula.format('abc')

    
    def read_atoms_collection(self, ):
        """
        """
        atoms = Atoms()
        # atoms
        coll_atom_kinds = [coll for coll in coll.children if 'atoms' in coll.name][0]
        coll_cell = [coll for coll in coll.children if 'cell' in coll.name][0]
        for obj in coll_atom_kinds.all_objects:
            print(obj.name)
            ind = obj.name.index('atom_kind_')
            ele = obj.name[ind + 10:].split('_')[0].split('.')[0]
            if len(obj.children) != 0:
                for vertex in obj.data.vertices:
                    location = obj.matrix_world @ vertex.co
                    atoms.append(Atom(ele, location))
            else:
                if not obj.parent:
                    location = obj.location
                    atoms.append(Atom(ele, location))
        # cell
        cell_vertexs = []
        if 'point_cell' in coll_cell.all_objects.keys():
            obj = coll_cell.all_objects['point_cell']
            for vertex in obj.data.vertices:
                location = obj.matrix_world @ vertex.co
                cell_vertexs.append(location)
        print(atoms)
        print('cell_vertexs: ', cell_vertexs)
        if cell_vertexs:
            cell = [cell_vertexs[4], cell_vertexs[2], cell_vertexs[1]]
            atoms.cell = cell
            atoms.pbc = [True, True, True]
            print('cell: ', cell)
        # self.atoms = atoms
        return atoms
    def draw(self, ):
        update = False
        for obj in bpy.context.selected_objects:
            if self.name == obj.name:
                update = True
        if update:
            self.update()
        else:
            self._draw()
    def _draw(self, ):
        obj = Blase(self)
        obj.draw()
        self.bobj = obj
    def update(self, ):
        self.read_atoms_collection()
        # self.atoms = atoms
        pass
        
    def delete(self, index = []):
        """
        delete atoms
        (1) atoms kinds, positions
        (2) obj_atom
        """
        # atoms.bobj.coll.children[1].all_objects[0].data.vertices[0].co
        for kind, mesh in self.bobj.meshes.items():
            bm = bmesh.new()
            bm.from_mesh(mesh)
            bm.verts.ensure_lookup_table()
            inds = [atom.index for atom in self if self.kinds[atom.index]==kind]
            verts = []
            for i in range(len(bm.verts)):
                if inds[i] in index:
                    bm.verts.remove(bm.verts[i])
            bm.to_mesh(mesh)
        self.__delitem__(index)
    def translate(self, displacement):
        """Translate atomic positions.

        The displacement argument can be a float an xyz vector or an
        nx3 array (where n is the number of atoms)."""

        self.arrays['positions'] += np.array(displacement)
        for kind, mesh in self.bobj.meshes.items():
            for vertice in mesh.vertices:
                vertice.co += mu.Vector(displacement)

    def choose_objects(action_type,
                   model_type):
        """
        """
        # For selected objects of all selected layers
        change_objects = []
        # Note all selected objects first.
        for atom in bpy.context.selected_objects:
            change_objects.append(atom)
        print(change_objects)

        # This is very important now: If there are dupliverts structures, note
        # only the parents and NOT the children! Otherwise the double work is
        # done or the system can even crash if objects are deleted. - The
        # chidlren are accessed anyways (see below).
        change_objects = []
        for atom in change_objects_all:
            if atom.parent != None:
                FLAG = False
                for atom2 in change_objects:
                    if atom2 == atom.parent:
                       FLAG = True
                if FLAG == False:
                    change_objects.append(atom)
            else:
                change_objects.append(atom)

        # And now, consider all objects, which are in the list 'change_objects'.
        for atom in change_objects:
            if len(atom.children) != 0:
                for atom_child in atom.children:
                    if atom_child.type in {'SURFACE', 'MESH', 'META'}:
                        modify_objects(action_type,
                                       atom_child,
                                       model_type)
            else:
                if atom.type in {'SURFACE', 'MESH', 'META'}:
                    modify_objects(action_type,
                                   atom,
                                   model_type)
