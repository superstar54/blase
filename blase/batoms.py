"""Definition of the Atoms class.

This module defines the atoms object in the blase package.

"""

from ase import Atom, Atoms
import os
import bpy
import mathutils as mu
import bmesh
from mathutils import Vector
from math import sqrt
from copy import copy
from blase.tools import get_atom_kinds, get_bond_kinds, get_bondpairs, get_cell_vertices, get_polyhedra_kinds, search_pbc
from blase.btools import draw_cell, draw_atoms, draw_atom_kind, draw_bonds, draw_polyhedras, draw_isosurface, bond_source, cylinder_mesh_from_instance, clean_default
import numpy as np
from ase.cell import Cell




subcollections = ['atoms', 'bonds', 'instancers', 'cell', 'polyhedras', 'isosurfaces', 'boundary']




class Batoms():
    """Atoms object.

    Parameters:

    
    """

    def __init__(self, atoms, name = None, coll = None, draw = False, model_type = '0',
                 scale = 1.0):
        if not hasattr(bpy.types.Collection, 'blase'):
            print('add blase')
            bpy.types.Collection.blase = bpy.props.PointerProperty(name="blase", type = BlaseSettings)
            # bpy.types.Collection.batoms = bpy.props.CollectionProperty(name="batoms", type = BlaseAtom)
        self.atoms = atoms
        self.name = name
        if not self.name:
            self.name = self.atoms.symbols.formula.format('abc')
        self.scene = bpy.context.scene
        self.coll = coll
        self.set_collection()
        # self.set_atoms()
        self.coll.blase.model_type = model_type
        self.scale = scale
        #
        
        if draw:
            self.draw()
        

    # def set_atoms(self):
    #     for atom in self.atoms:
    #         batom = self.coll.batoms.add()
    #         batom.symbol = atom.symbol
    #         batom.position = atom.position
    def set_collection(self):
        if not self.coll:
            for coll in bpy.data.collections:
                if self.name == coll.name:
                    raise Exception("Failed, the name %s already in use!"%self.name)
            self.coll = bpy.data.collections.new(self.name)
            self.coll.blase.is_blase = True
            self.scene.collection.children.link(self.coll)
            for sub_name in subcollections:
                subcoll = bpy.data.collections.new('%s_%s'%(self.name, sub_name))
                self.coll.children.link(subcoll)
        
    def draw_cell(self):
        # bond_kinds
        # if not self.bond_list and self.bond_cutoff:
        cell_vertices = get_cell_vertices(self.atoms.cell)
        for obj in self.coll.children['%s_cell'%self.name].all_objects:
            bpy.data.objects.remove(obj)
        draw_cell(self.coll.children['%s_cell'%self.name], cell_vertices)
    def draw_atoms(self, scale = 1.0, props = {}):
        # atom_kinds
        for obj in self.coll.children['%s_atoms'%self.name].all_objects:
            bpy.data.objects.remove(obj)
        for obj in self.coll.children['%s_instancers'%self.name].all_objects:
            bpy.data.objects.remove(obj)
        self.scale = scale
        self.atom_kinds = get_atom_kinds(self.atoms, scale = self.scale, props = props)
        self.nkinds = len(self.atom_kinds)
        draw_atoms(self.coll.children['%s_atoms'%self.name], self.atom_kinds)
        self.draw_atoms_boundary()
    def draw_atoms_boundary(self, cutoff = None, props = {}):
        # atom_kinds
        if not cutoff:
            cutoff = self.coll.blase.boundary
        self.bd_atoms, self.bd_index = search_pbc(self.atoms, cutoff)
        for obj in self.coll.children['%s_boundary'%self.name].all_objects:
            bpy.data.objects.remove(obj)
        self.bd_atom_kinds = get_atom_kinds(self.bd_atoms, scale = self.scale, props = props)
        self.nkinds = len(self.bd_atom_kinds)
        draw_atoms(self.coll.children['%s_boundary'%self.name], self.bd_atom_kinds)
        atoms = self.atoms + self.bd_atoms
        # if self.coll.blase.model_type in ['1', '2', '3']:
            # self.draw_bonds(atoms = atoms)
    def draw_bonds(self, cutoff = 1.0, atoms = None):
        # bond_kinds
        # if not self.bond_list and self.bond_cutoff:
        for obj in self.coll.children['%s_bonds'%self.name].all_objects:
            bpy.data.objects.remove(obj)
        if not atoms:
            atoms = self.atoms
        self.bond_list = get_bondpairs(atoms)
        self.bond_kinds = get_bond_kinds(atoms, self.bond_list)
        draw_bonds(self.coll.children['%s_bonds'%self.name], self.bond_kinds)
    def draw_isosurface(self):
        if self.isosurface:
            volume = self.isosurface[0]
            icolor = 0
            if len(self.isosurface) == 1:
                draw_isosurface(self.coll.children['%s_isosurface'%self.name], volume, level=None, icolor = icolor)
            for level in self.isosurface[1:]:
                draw_isosurface(self.coll.children['%s_isosurface'%self.name], volume, level=level, icolor = icolor)
                icolor += 1
    def draw(self, model_type = None):
        if not model_type:
            model_type = self.coll.blase.model_type
        if model_type == '0':
            self.draw_cell()
            self.draw_atoms()
        elif model_type == '1':
            # view(images)
            self.draw_cell()
            self.draw_atoms(scale = 0.6)
            self.draw_bonds()
        elif model_type == '2':
            self.draw_cell()
            self.draw_atoms(scale = 0.6)
            self.draw_bonds()
            self.draw_polyhedras()
        elif model_type == '3':
            self.draw_cell()
            self.draw_atoms(scale = 0.001)
            self.draw_bonds()
        # self.draw_atoms()
        # self.draw_cell()
        # self.draw_bonds()
    def replace(self, index = [], element = None):
        """
        replace atoms
        """
        from blase.tools import default_atom_kind
        # delete old verts
        self.delete_verts(index)
        # if kind exists, merger, otherwise build a new kind and add.
        name = 'atom_kind_{0}'.format(element)
        positions = self.atoms.positions[index]
        if name in self.coll.children['%s_atoms'%self.name].objects:
            obj = self.coll.children['%s_atoms'%self.name].objects[name]
            bm = bmesh.new()
            bm.from_mesh(obj.data)
            bm.verts.ensure_lookup_table()
            verts = []
            for pos in positions:
                bm.verts.new(pos)
            bm.to_mesh(obj.data)
        else:
            atom_kind = default_atom_kind(element, positions)
            draw_atom_kind(element, self.coll.children['%s_atoms'%self.name], atom_kind)
        # ase atoms
        self.atoms.symbols[index] = element
    def delete_verts(self, index = []):
        """
        delete atoms
        (1) atoms kinds, positions
        (2) obj_atom
        """
        # atoms.bobj.coll.children[1].all_objects[0].data.vertices[0].co
        for kind, obj in self.coll.children['%s_atoms'%self.name].objects.items():
            bm = bmesh.new()
            bm.from_mesh(obj.data)
            bm.verts.ensure_lookup_table()
            inds = [atom.index for atom in self.atoms if atom.symbol==kind.split('_')[-1]]
            verts_select = [bm.verts[i] for i in range(len(bm.verts)) if inds[i] in index] 
            bmesh.ops.delete(bm, geom=verts_select, context='VERTS')
            if len(bm.verts) == 0:
                bpy.data.objects.remove(obj)
            else:
                bm.to_mesh(obj.data)
    def delete(self, index = []):
        self.delete_verts(index)
        del self.atoms[index]
    def translate(self, displacement):
        """Translate atomic positions.

        The displacement argument can be a float an xyz vector or an
        nx3 array (where n is the number of atoms)."""

        self.atoms.arrays['positions'] += np.array(displacement)
        for kind, obj in self.coll.children['%s_atoms'%self.name].objects.items():
            for vertice in obj.data.vertices:
                vertice.co += mu.Vector(displacement)
    def __add__(self, other):
        self += other
        return self
    def __iadd__(self, other):
        self.extend(other)
        return self
    def extend(self, other):
        """Extend atoms object by appending atoms from *other*."""
        self.atoms = self.atoms + other.atoms
        for kind, obj in other.coll.children['%s_atoms'%other.name].objects.items():
            if kind in self.coll.children['%s_atoms'%self.name].objects:
                objs = [obj, self.coll.children['%s_atoms'%self.name].objects[kind]]
                bpy.ops.object.join(objs)
            else:
                self.coll.children['%s_atoms'%self.name].objects.link(obj)
        bpy.data.collections.remove(other.coll)
    def __imul__(self, m):
        """In-place repeat of atoms."""
        self.atoms *= m
        if 'kinds' in self.atoms.info:
            del self.atoms.info['kinds']
        self.update()
        return self
    def repeat(self, rep):
        """Create new repeated atoms object."""
        self.__imul__(rep)
    def __mul__(self, rep):
        self.repeat(rep)
        return self
    def copy(self, name = None):
        """Return a copy."""
        if not name:
            name = self.name + 'copy'
        batoms = self.__class__(self.atoms, name = name, draw = True)
        return batoms
    def update(self):
        """
        """
        self.draw()

    def set_cell(self, cell, scale_atoms=False):
        """Set unit cell vectors.

        Parameters:

        Examples:

        """

        # Override pbcs if and only if given a Cell object:
        self.atoms.set_cell(cell, scale_atoms = scale_atoms)
        self.update()
    
    def draw_constraints(self):
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
        coll_highlight = bpy.data.collections.new('highlight')
        self.coll.children.link(coll_highlight)
        # build materials
        material = bpy.data.materials.new('highlight')
        material.name = 'highlight'
        material.diffuse_color = color + (transmit,)
        # material.alpha_threshold = 0.2
        material.blend_method = 'BLEND'
        material.use_nodes = True
        principled_node = material.node_tree.nodes['Principled BSDF']
        principled_node.inputs['Alpha'].default_value = transmit
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
    def load_frames(self):
        # render settings
        for i in range(0, self.nimages):
            atom_kinds = get_atom_kinds(self.images[i])
            # bond_kinds = get_bond_kinds(self.images[i], self.bond_list, self.nbins)
            for kind, datas in atom_kinds.items():
                obj_atom = bpy.data.objects['atom_kind_{0}'.format(kind)]
                nverts = len(obj_atom.data.vertices)
                for j in range(nverts):
                    obj_atom.data.vertices[j].co = datas['positions'][j]
                    obj_atom.data.vertices[j].keyframe_insert('co', frame=i + 1)
        self.scene.frame_start = 1
        self.scene.frame_end = self.nimages
    

