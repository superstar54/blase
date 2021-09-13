"""Definition of the Batom class.

This module defines the Batom object in the blase package.

"""

from ase import Atom, Atoms
import bpy
import bmesh
from mathutils import Vector
from copy import copy
from blase.data import material_styles_dict
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
                label,
                species,
                positions,
                element = None,
                scale = 1.0, 
                kind_props = {},
                color_style = 'JMOL',
                material_style = 'blase',
                bsdf_inputs = None,
                draw = False, 
                 ):
        #
        self.label = label
        self.species = species
        if not element:
            self.element = species.split('_')[0]
        else:
            self.element = element
        self.scene = bpy.context.scene
        self.name = 'atom_%s_%s'%(self.label, self.species)
        self.kind_props = kind_props
        self.color_style = color_style
        self.material_style = material_style
        self.bsdf_inputs = bsdf_inputs
        self.species_data = get_atom_kind(self.element, scale = scale, props = self.kind_props, color_style = self.color_style)
        self.set_material()
        self.set_instancer()
        self.set_object(positions)
        self.bond_data = {}
        self.polyhedra_data = {}
        if draw:
            self.draw_atom()
    def set_material(self):
        if self.name not in bpy.data.materials:
            if not self.bsdf_inputs:
                bsdf_inputs = material_styles_dict[self.material_style]
            material = bpy.data.materials.new(self.name)
            material.diffuse_color = np.append(self.species_data['color'], self.species_data['transmit'])
            material.metallic = bsdf_inputs['Metallic']
            material.roughness = bsdf_inputs['Roughness']
            material.blend_method = 'BLEND'
            material.use_nodes = True
            principled_node = material.node_tree.nodes['Principled BSDF']
            principled_node.inputs['Base Color'].default_value = np.append(self.species_data['color'], self.species_data['transmit'])
            principled_node.inputs['Alpha'].default_value = self.species_data['transmit']
            for key, value in bsdf_inputs.items():
                principled_node.inputs[key].default_value = value
        else:
            material = bpy.data.materials[self.name]
        self.material = material
    def set_instancer(self):
        name = 'instancer_atom_{0}_{1}'.format(self.label, self.species)
        if name not in bpy.data.objects:
            bpy.ops.mesh.primitive_uv_sphere_add(radius = self.species_data['radius']) #, segments=32, ring_count=16)
            sphere = bpy.context.view_layer.objects.active
            if isinstance(self.species_data['scale'], float):
                self.species_data['scale'] = [self.species_data['scale']]*3
            sphere.scale = self.species_data['scale']
            sphere.name = 'instancer_atom_{0}_{1}'.format(self.label, self.species)
            sphere.data.materials.append(self.material)
            bpy.ops.object.shade_smooth()
            sphere.hide_set(True)
        else:
            sphere = bpy.data.objects[name]
        self.instancer = sphere
    def set_object(self, positions):
        """
        build child object and add it to main objects.
        """
        if self.name not in bpy.data.objects:
            mesh = bpy.data.meshes.new(self.name)
            obj_atom = bpy.data.objects.new(self.name, mesh)
            obj_atom.data.from_pydata(positions, [], [])
            obj_atom.is_batom = True
        elif hasattr(bpy.data.objects[self.name], 'batom'):
            obj_atom = bpy.data.objects[self.name]
        else:
            raise Exception("Failed, the name %s already in use and is not blase object!"%self.name)
        obj_atom.species = self.species
        obj_atom.element = self.element
        self.instancer.parent = obj_atom
        obj_atom.instance_type = 'VERTS'
        bpy.data.collections['Collection'].objects.link(obj_atom)
        self.batom = obj_atom
    @property
    def scale(self):
        return self.get_scale()
    @scale.setter
    def scale(self, scale):
        self.set_scale(scale)
    def get_scale(self):
        return self.batom.scale
    def set_scale(self, scale):
        if isinstance(scale, float) or isinstance(scale, int):
            scale = [scale]*3
        self.instancer.scale = scale
    @property
    def positions(self):
        return self.get_positions()
    @positions.setter
    def positions(self, positions):
        self.set_positions(positions)
    def get_positions(self):
        """
        Get array of positions.
        """
        return np.array([self.batom.data.vertices[i].co for i in range(len(self))])
    def set_positions(self, positions):
        """
        Set positions
        """
        natoms = len(self)
        if len(positions) != natoms:
            raise ValueError('positions has wrong shape %s != %s.' %
                                (len(positions), natoms))
        for i in range(natoms):
            self.batom.data.vertices[i].co = positions[i]
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
        obj = self.batom
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
    def __delitem__(self, index):
        """
        """
        self.delete(index)
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
    
    def __len__(self):
        return len(self.batom.data.vertices)
    
    def __getitem__(self, index):
        """Return a subset of the Batom.

        i -- int, describing which atom to return.
        """

        if isinstance(index, int):
            natoms = len(self)
            if index < -natoms or index >= natoms:
                raise IndexError('Index out of range.')
            return self.batom.data.vertices[index].co
    def __setitem__(self, i, value):
        """Return a subset of the Batom.

        i -- int, describing which atom to return.
        """
        if isinstance(i, int):
            natoms = len(self)
            if i < -natoms or i >= natoms:
                raise IndexError('Index out of range.')
            self.batom.data.vertices[i].co = value

    def repeat(self, m, cell):
        """
        In-place repeat of atoms.
        """
        if isinstance(m, int):
            m = (m, m, m)
        for x, vec in zip(m, cell):
            if x != 1 and not vec.any():
                raise ValueError('Cannot repeat along undefined lattice '
                                 'vector')
        M = np.product(m)
        n = len(self)
        
        self.positions = self.positions
        positions = np.tile(self.positions, (M,) + (1,) * (len(self.positions.shape) - 1))
        i0 = 0
        for m0 in range(m[0]):
            for m1 in range(m[1]):
                for m2 in range(m[2]):
                    i1 = i0 + n
                    positions[i0:i1] += np.dot((m0, m1, m2), cell)
                    i0 = i1
        self.add_vertices(positions[n:])
    def copy(self, label, species):
        """
        Return a copy.

        name: str
            The name of the copy.

        For example, copy H species:
        
        >>> h_new = h.copy(name = 'h_new')

        """
        name = 'atom_%s_%s'%(label, species)
        obj = bpy.data.objects.new(name, self.batom.data)
        name = 'instancer_atom_{0}_{1}'.format(self.label, self.species)
        sphere = bpy.data.objects.new(name, self.instancer.data)
        mat = self.material.copy()
        mat.name = name
        batom = self.__class__(label, species, self.positions)
        return batom
    def extend(self, other):
        """
        Extend batom object by appending batom from *other*.
        
        >>> from blase.batoms import Batom
        >>> h = Batom('h2o', 'H', [[0, 0, 0], [2, 0, 0]])
        >>> o = Batom('h2o', 'O', [[0, 0, 0]])
        >>> o.extend(h)
        """
        bpy.ops.object.select_all(action='DESELECT')
        self.batom.select_set(True)
        other.batom.select_set(True)
        bpy.context.view_layer.objects.active = self.batom
        bpy.ops.object.join()
    def add_vertices(self, positions):
        """
        """
        bm = bmesh.new()
        bm.from_mesh(self.batom.data)
        bm.verts.ensure_lookup_table()
        verts = []
        for pos in positions:
            bm.verts.new(pos)
        bm.to_mesh(self.batom.data)
    def translate(self, displacement):
        """Translate atomic positions.

        The displacement argument is an xyz vector.

        For example, move H species molecule by a vector [0, 0, 5]

        >>> h.translate([0, 0, 5])
        """
        bpy.ops.object.select_all(action='DESELECT')
        self.batom.select_set(True)
        bpy.ops.transform.translate(value=displacement)
    def rotate(self, angle, axis = 'Z', orient_type = 'GLOBAL'):
        """Rotate atomic based on a axis and an angle.

        Parameters:

        angle: float
            Angle that the atoms is rotated around the axis.
        axis: str
            'X', 'Y' or 'Z'.

        For example, rotate h2o molecule 90 degree around 'Z' axis:
        
        >>> h.rotate(90, 'Z')

        """
        bpy.ops.object.select_all(action='DESELECT')
        self.batom.select_set(True)
        bpy.ops.transform.rotate(value=angle, orient_axis=axis.upper(), orient_type = orient_type)
    
    
    