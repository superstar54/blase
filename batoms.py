"""Definition of the Batoms class.

This module defines the batoms object in the blase package.

"""

from ase import Atom, Atoms
import bpy
import bmesh
from mathutils import Vector
from copy import copy
from blase.tools import get_atom_kinds, get_bond_kinds, get_bondpairs, get_cell_vertices, get_polyhedra_kinds, search_pbc, get_bbox
from blase.bdraw import draw_cell, draw_atoms, draw_atom_kind, draw_bonds, draw_polyhedras, draw_isosurface, bond_source, cylinder_mesh_from_instance, clean_default
import numpy as np
import time

subcollections = ['atoms', 'bonds', 'instancers', 'cell', 'polyhedras', 'isosurface', 'boundary']


class Batoms():
    """Batoms Class

    In Blender, the atoms objects and collections are organised in the following way, 
    take water molecule as a example:

    * h2o                            # main collection

      * h2o_atoms                   # atoms collection
    
        * atom_h2o_H                  # atoms object
    
        * atom_h2o_O                  # atoms object

      * h2o_bonds                    # bonds collection

        * bond_h2o_H                 --bond object
    
        * bond_h2o_O                 --bond object
    
      * h2o_cell                     --cell collection
    
      * h2o_instancers               --instancer collection
    
        * sphere_h2o_H               --sphere object
    
        * sphere_h2o_O               --sphere object
    
    Then, a Batoms object is linked to this main collection in Blender. 
    

    Parameters:

    atoms: ase.atoms.Atoms object
        The atomic structure built using ASE
    coll: bpy.types.Collection
        The main Blender collection. 
        This collection has a ``blase`` property, which has five parameters:
            * is_blase
            * pbc
            * unit cell
            * boundary
            * model_type
    model_type: str
        enum in '0', '1', '2', '3'
    boundary:  list 
        search atoms at the boundary
    add_bonds: dict
        add bonds not in the default
    
    Examples:

    >>> from ase.build import molecule
    >>> from blase.batoms import Batoms
    >>> h2o = molecule('H2O')
    >>> h2o = Batoms(h2o)
    >>> h2o.draw()
    """
    

    def __init__(self, 
                 coll = None,
                 atoms = None, 
                 name = None,
                 model_type = '0', 
                 scale = 1.0, 
                 boundary = [0.0, 0.0, 0.0],
                 bond_list = None,
                 add_bonds = {}, 
                 remove_bonds = {},
                 hydrogen_bond = None,
                 polyhedra_dict = {},
                 show_unit_cell = True,
                 isosurface = [],
                 kind_props = {},
                 color = 'JMOL',
                 movie = False,
                 draw = False, 
                 ):
        #
        self.scene = bpy.context.scene
        self.scale = scale
        self.bond_list = bond_list
        self.add_bonds = add_bonds
        self.remove_bonds = remove_bonds
        self.hydrogen_bond = hydrogen_bond
        self.polyhedra_dict = polyhedra_dict
        self.show_unit_cell = show_unit_cell
        self.isosurface = isosurface
        self.kind_props = kind_props
        self.name = name
        self.color = color
        if atoms:
            if not isinstance(atoms, list):
                atoms = [atoms]
            self.images = atoms
            # batoms save the real objects
            self.batoms = atoms[0]
            # atms save all the atoms to be drawed. including the boundary atoms.
            self.atoms = self.batoms.copy()
            if not self.name:
                self.name = self.batoms.symbols.formula.format('abc')
            # add new collection
            self.set_collection()
            # save atoms to blase.atoms
            # self.atoms2batoms()
            self.coll.blase.model_type = model_type
            self.coll.blase.boundary = boundary
            self.coll.blase.pbc = self.npbool2bool(self.batoms.pbc)
            self.coll.blase.cell = self.batoms.cell[:].flatten()
        elif coll:
            print('Build from collection')
            self.coll = coll
            self.coll2atoms()
            # self.batoms = self.batoms2atoms()
            self.atoms = self.batoms.copy()
        else:
            raise Exception("Failed, either atoms or coll should be provided!"%self.name)
        if draw:
            self.draw()
        if movie:
            self.load_frames()
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
    def atoms2batoms(self, atoms = None):
        """
        """
        tstart = time.time()
        if not atoms:
            atoms = self.atoms
        for atom in atoms:
            batom = self.coll.batoms.add()
            batom.symbol = atom.symbol
            batom.position = atom.position
        print('set atoms2batoms: {0:10.2f} s'.format(time.time() - tstart))
    def coll2atoms(self):
        """
        build ASE Atoms from a blase's main collection
        """
        atoms = Atoms()
        coll = self.coll
        # atoms
        for obj in coll.children['%s_atoms'%coll.name].all_objects:
            ele = obj.name.split('_')[2]
            if len(obj.children) != 0:
                for vertex in obj.data.vertices:
                    location = obj.matrix_world @ vertex.co
                    atoms.append(Atom(symbol = ele, position = location))
            else:
                if not obj.parent:
                    location = obj.location
                    atoms.append(Atom(ele, location))
        # cell
        coll_cell = coll.children['%s_cell'%coll.name]
        cell_vertexs = []
        if 'point_cell' in coll_cell.all_objects.keys():
            obj = coll_cell.all_objects['point_cell']
            for vertex in obj.data.vertices:
                location = obj.matrix_world @ vertex.co
                cell_vertexs.append(location)
        if cell_vertexs:
            cell = [cell_vertexs[4], cell_vertexs[2], cell_vertexs[1]]
            atoms.cell = cell
            atoms.pbc = coll.blase.pbc
        self.atoms = atoms
        self.batoms = atoms
    def batoms2atoms(self, coll = None):
        """
        """
        from ase import Atom, Atoms
        atoms = Atoms()
        tstart = time.time()
        if not coll:
            coll = self.coll
        for batom in coll.batoms:
            atom = Atom(symbol = batom.symbol, 
                        position = batom.position)
            atoms.append(atom)
        atoms.pbc = coll.blase.pbc
        atoms.cell = np.array(coll.blase.cell).reshape(3, 3)
        print('set batoms2atoms: {0:10.2f} s'.format(time.time() - tstart))
        return atoms
    def set_collection(self):
        """
        build main collection and its child collections.
        """
        for coll in bpy.data.collections:
            if self.name == coll.name:
                raise Exception("Failed, the name %s already in use!"%self.name)
        self.coll = bpy.data.collections.new(self.name)
        self.coll.blase.is_blase = True
        self.scene.collection.children.link(self.coll)
        for sub_name in subcollections:
            subcoll = bpy.data.collections.new('%s_%s'%(self.name, sub_name))
            self.coll.children.link(subcoll)
    def search_boundary(self, ):
        """
        Search pbc
        """
        boundary =  np.array(self.coll.blase.boundary[:]) > 0.0
        if self.atoms.pbc.any() and boundary.any():
            bd_atoms, bd_index = search_pbc(self.atoms, self.coll.blase.boundary)
            self.atoms = self.batoms + bd_atoms
            self.atoms.info['species'] = self.atoms.info['species'] + bd_atoms.info['species']
        else:
            self.atoms = self.batoms
    def draw_cell(self):
        """
        Draw unit cell
        """
        print('--------------Draw cell--------------')
        cell_vertices = get_cell_vertices(self.coll.blase.cell)
        for obj in self.coll.children['%s_cell'%self.name].all_objects:
            bpy.data.objects.remove(obj)
        if self.show_unit_cell:
            draw_cell(self.coll.children['%s_cell'%self.name], cell_vertices)
    def draw_atoms(self, scale = None, kind_props = {}):
        """
        Draw atoms.

        Parameters:

        scale: float
            scale for the sphere
        kind_props: dict
            Set user defined properties for species
        """
        print('--------------Draw atoms--------------')
        self.clean_atoms()
        if scale:
            self.scale = scale
        if not kind_props:
            kind_props = self.kind_props
        self.search_boundary()
        self.atom_kinds = get_atom_kinds(self.atoms, scale = self.scale, props = kind_props, color = self.color)
        draw_atoms(self.coll.children['%s_atoms'%self.name], self.atom_kinds)
    def draw_bonds(self, cutoff = 1.0):
        """
        Draw bonds.
        Parameters:

        cutoff: float
            cutoff used to build bond pairs.
        """
        print('--------------Draw bonds--------------')
        self.clean_bonds()
        # print(self.atoms.info)
        if not self.bond_list:
            self.bond_list = get_bondpairs(self.atoms, add_bonds = self.add_bonds, remove_bonds=self.remove_bonds)
        if self.hydrogen_bond:
            self.hydrogen_bond_list = get_bondpairs(self.atoms, cutoff = {('O', 'H'): self.hydrogen_bond})
            self.bond_kinds = get_bond_kinds(self.atoms, self.atom_kinds, self.hydrogen_bond_list, color = self.color)
            draw_bonds(self.coll.children['%s_bonds'%self.name], self.bond_kinds)
        self.bond_kinds = get_bond_kinds(self.atoms, self.atom_kinds, self.bond_list, color = self.color)
        draw_bonds(self.coll.children['%s_bonds'%self.name], self.bond_kinds)
        
    def draw_polyhedras(self, cutoff = 1.0):
        """
        Draw bonds.
        Parameters:

        cutoff: float
            cutoff used to build bond pairs.
        """
        print('--------------Draw bonds--------------')
        self.clean_polyhedras()
        polyhedra_kinds = get_polyhedra_kinds(self.atoms, self.atom_kinds, self.bond_list, polyhedra_dict = self.polyhedra_dict)
        draw_polyhedras(self.coll.children['%s_polyhedras'%self.name], polyhedra_kinds)
    def draw_isosurface(self, isosurface = []):
        """
        Draw bonds.

        Parameters:

        isosurface: list
            isosurface data.
        """
        if not isosurface:
            isosurface = self.isosurface
        volume = self.isosurface[0]
        icolor = 0
        if len(self.isosurface) == 1:
            draw_isosurface(self.coll.children['%s_isosurface'%self.name], volume, cell = self.atoms.cell, level=None, icolor = icolor)
        for level in self.isosurface[1:]:
            draw_isosurface(self.coll.children['%s_isosurface'%self.name], volume, cell = self.atoms.cell, level=level, icolor = icolor)
            icolor += 1
    def draw(self, model_type = None):
        """
        Draw atoms, bonds, polyhedra.

        Parameters:

        model_type: str

        """
        if not model_type:
            model_type = self.coll.blase.model_type
        else:
            self.coll.blase.model_type = model_type
        self.draw_cell()
        if model_type == '0':
            self.scale = 1.0
            self.draw_atoms()
            self.clean_bonds()
        elif model_type == '1':
            # view(images)
            self.scale = 0.4
            self.draw_atoms()
            self.draw_bonds()
        elif model_type == '2':
            self.scale = 0.4
            self.draw_atoms()
            self.draw_bonds()
            self.draw_polyhedras()
        elif model_type == '3':
            self.scale = 0.01
            self.draw_atoms()
            # self.clean_atoms
            self.draw_bonds()
        if self.isosurface:
            self.draw_isosurface()
    def clean_bonds(self, ):
        """
        remove all bond object in the bond collection
        """
        for obj in self.coll.children['%s_bonds'%self.name].all_objects:
            bpy.data.objects.remove(obj)
    def clean_atoms(self, ):
        """
        remove all atom object in the atom collection
        """
        for obj in self.coll.children['%s_atoms'%self.name].all_objects:
            bpy.data.objects.remove(obj)
        for obj in self.coll.children['%s_instancers'%self.name].all_objects:
            bpy.data.objects.remove(obj)
    def clean_polyhedras(self, ):
        """
        remove all polyhedra object in the polyhedra collection
        """
        for obj in self.coll.children['%s_polyhedras'%self.name].all_objects:
            bpy.data.objects.remove(obj)
    def replace(self, index = [], element = None):
        """
        replace atoms.

        Parameters:

        index: list
            index of atoms will be replaced.

        element: str
            atoms will be changed to this element.
        
        >>> h2o.replace([1, 2], 'S')

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
            atom_kind = default_atom_kind(element, positions, color=self.color)
            draw_atom_kind(element, self.coll.children['%s_atoms'%self.name], atom_kind)
        # ase atoms
        self.atoms.symbols[index] = element
    def delete_verts(self, index = []):
        """
        delete verts
        
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
        """
        delete atoms.

        index: list
            index of atoms to be delete
        
        For example, delete the second atom in h2o molecule. 
        Please note that index start from 0.

        >>> h2o.delete([1])

        """
        self.delete_verts(index)
        del self.atoms[index]
    def translate(self, displacement):
        """Translate atomic positions.

        The displacement argument is an xyz vector.

        For example, move h2o molecule by a vector [0, 0, 5]

        >>> h2o.translate([0, 0, 5])

        """
        self.update_collection()
        self.atoms.translate(displacement)
        bpy.ops.object.select_all(action='DESELECT')
        for coll in self.coll.children:
            for obj in coll.objects:
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
        self.update_collection()
        self.atoms.rotate(angle, axis.lower())
        bpy.ops.object.select_all(action='DESELECT')
        for coll in self.coll.children:
            for obj in coll.objects:
                obj.select_set(True)
        bpy.ops.transform.rotate(value=angle, orient_axis=axis.upper(), orient_type = orient_type)
        
    def __add__(self, other):
        self += other
        return self
    def __iadd__(self, other):
        self.extend(other)
        return self
    def extend(self, other):
        """
        Extend atoms object by appending atoms from *other*.
        
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
        self.atoms = self.atoms + other.atoms
        for kind, obj in other.coll.children['%s_atoms'%other.name].objects.items():
            if kind in self.coll.children['%s_atoms'%self.name].objects:
                objs = [obj, self.coll.children['%s_atoms'%self.name].objects[kind]]
                bpy.ops.object.join(objs)
            else:
                self.coll.children['%s_atoms'%self.name].objects.link(obj)
        bpy.data.collections.remove(other.coll)
    def __imul__(self, m):
        """
        In-place repeat of atoms.
        """
        print(m)
        self.atoms *= m
        if 'species' in self.atoms.info:
            del self.atoms.info['species']
        self.update()
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
    def copy(self, name = None):
        """
        Return a copy.

        name: str
            The name of the copy.

        For example, copy h2o molecule:
        
        >>> h2o_new = h2o.copy(name = 'h2o_new')

        """
        if not name:
            name = self.name + 'copy'
        batoms = self.__class__(atoms = self.atoms, name = name, model_type = self.coll.blase.model_type, draw = True)
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
        self.draw()
    def update_collection(self):
        """
        """
        self.coll = bpy.data.collections[self.name]
        self.coll2atoms()
        return 0


    def set_cell(self, cell, scale_atoms=False):
        """Set unit cell vectors.

        Parameters:

        Examples:

        """

        # Override pbcs if and only if given a Cell object:
        self.atoms.set_cell(cell, scale_atoms = scale_atoms)
        self.coll.blase.cell = self.atoms.cell[:].flatten()
        self.update()
    
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
    
    def render(self, **kwargs):
        """
        Render the atoms, and save to a png image.

        Support all parameters for Class Blase

        >>> h2o.render(resolution_x = 1000, output_image = 'h2o.png')
        
        """
        from blase.bio import Blase
        print('--------------Render--------------')
        print('Rendering atoms')
        if 'bbox' not in kwargs:
            bbox = get_bbox(bbox = None, atoms = self.atoms)
            kwargs['bbox'] = bbox
        if 'output_image' not in kwargs:
            kwargs['output_image'] = '%s.png'%self.name
        # print(kwargs)
        bobj = Blase(**kwargs)
        for function in bobj.functions:
            name, paras = function
            getattr(bobj, name)(**paras)
        # bobj.load_frames()
        bobj.render()
    def set_model_type(self, model_type = 1) -> None:
        """
        set model type
        """
        self.update_collection()
        self.coll.blase.model_type = str(model_type)
        self.draw()
