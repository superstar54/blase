# ###################################################################
"""Python module for drawing and rendering ase atoms objects using blender.

Part of code based on ase.io.pov, 
https://wiki.fysik.dtu.dk/ase/dev/_modules/ase/io/pov.html

"""

# ###################################################################

import bpy
from mathutils import Vector, Matrix
import os
import numpy as np
from ase import Atoms, Atom
from ase.io.utils import PlottingVariables, cell_to_lines
from ase.constraints import FixAtoms
from ase.utils import basestring
from ase.io import read, write
from ase.build import molecule, bulk
from ase.data import covalent_radii, atomic_numbers
from ase.data.colors import jmol_colors
from math import pi, sqrt, radians, acos, atan2
from blase.tools import get_atom_kinds, get_bond_kinds, get_bondpairs, get_polyhedra_kinds
from blase import tools
from blase.btools import draw_cell, draw_atoms, draw_bonds, draw_bonds_2, draw_polyhedras, draw_isosurface, bond_source, cylinder_mesh_from_instance, clean_default
from blase.connectivity import ConnectivityList
from blase.boundary import Boundary
import time



class Blase():
    """
    Blase object for drawing and rendering ase atoms.

    Examples:

    bobj = Blase(atoms, **kwargs0)
    bobj.draw_cell()
    bobj.draw_atoms()
    bobj.draw_bonds()
    bobj.render(filename)

    """
    default_settings = {
        
        'display': False,  # display while rendering
        'transparent': True,  # transparent background
        'resolution_x': None,  # 
        'resolution_y': None,  # 
        'camera': True,
        'camera_loc': None,  # x, y is the image plane, z is *out* of the screen
        'camera_type': 'ORTHO',  #  ['PERSP', 'ORTHO']
        'ortho_scale': None, #
        'camera_lens': 10,  #
        'fstop': 0.5,
        'camera_target': None, #
        'world': False,
        'light': True,
        'light_loc': [0, 0, 200],
        'light_type': 'SUN', # 'POINT', 'SUN', 'SPOT', 'AREA'
        'point_lights': [],  # 
        'light_strength': 1.0,
        'background': 'White',  # color
        'textures': None,  # length of atoms list of texture names
        'engine': 'BLENDER_EEVEE', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH'
        'transmits': None,  # transmittance of the atoms
        'show_unit_cell': 'default',
        'celllinewidth': 0.025,  # radius of the cylinders representing the cell
        'bbox': None,
        'bondlinewidth': 0.10,  # radius of the cylinders representing bonds
        'balltypes': None,
        'radii': None, 
        'colors': None,
        'make_real': False,
        'kind_props': None,
        'bond_cutoff': None,  # 
        'bond_list': {},  # [[atom1, atom2], ... ] pairs of bonding atoms
        'polyhedra_dict': {},
        'search_pbc_atoms': False, #{'bonds_dict': {}, 'molecule_list': {}},
        'search_molecule': False, #{'search_list': None},
        'boundary_list': [],
        'isosurface':None,
        'resolution_x': 1000,
        'cube': None,
        'highlight': None, # highlight atoms
        'functions': [],
        'run_render': True,
        'animation': False,
        'save_to_blend': False,
        'queue': None,
        'gpu': True,
        'num_samples': 128,
        'build_collection': True,
        }  
    #
    

    def __init__(self, images, outfile = 'bout', name = None, rotations=None, scale=1, debug = False,
                              **parameters):
        for k, v in self.default_settings.items():
            setattr(self, k, parameters.pop(k, v))
        #
        self.debug = debug
        if not isinstance(images, list):
            images = [images]
        self.nimages = len(images)
        self.outfile = outfile
        if rotations:
            for rotation in rotations:
                for i in range(self.nimages):
                    images[i].rotate(rotation[0], rotation[1], rotate_cell = True)
        #
        self.images = images
        self.atoms = images[0]
        if self.search_pbc_atoms:
            print(self.search_pbc_atoms)
            cl = ConnectivityList(self.atoms, cutoffs = self.bond_cutoff, **self.search_pbc_atoms)
            self.atoms = cl.build()
            # view(self.atoms)
            print(self.atoms)
        if self.boundary_list:
            bd = Boundary(self.atoms, self.boundary_list)
            self.atoms = bd.build()
        self.name = name
        if not self.name:
            self.name = self.atoms.symbols.formula.format('abc')
        self.natoms = len(self.atoms)
        self.positions = self.atoms.positions*scale
        self.numbers = self.atoms.get_atomic_numbers()
        self.symbols = self.atoms.get_chemical_symbols()
        self.cell = self.atoms.get_cell()*scale
        self.celllinewidth = self.celllinewidth*scale
        self.bondlinewidth = self.bondlinewidth*scale
        #------------------------------------------------------------
        # atom_kinds
        self.atom_kinds = get_atom_kinds(self.atoms, self.kind_props)
        self.nkinds = len(self.atom_kinds)
        if self.debug:
            print(self.atom_kinds)
        #
        if isinstance(self.colors, dict):
            for kind in self.colors:
                self.atom_kinds[kind]['color'] = self.colors[kind]
        if isinstance(self.balltypes, str):
            for kind in self.atom_kinds:
                self.atom_kinds[kind]['balltype'] = self.balltypes
        elif isinstance(self.balltypes, dict):
            for kind in self.balltypes:
                self.atom_kinds[kind]['balltype'] = self.balltypes[kind]
        #
        if isinstance(self.transmits, float):
            for kind in self.atom_kinds:
                self.atom_kinds[kind]['transmit'] = self.transmits
        elif isinstance(self.transmits, dict):
            for kind in self.transmits:
                self.atom_kinds[kind]['transmit'] = self.transmits[kind]
        if self.radii is None:
            self.radii = scale
        elif isinstance(self.radii, float):
            for kind in self.atom_kinds:
                self.atom_kinds[kind]['radius'] *= self.radii*scale
        elif isinstance(self.radii, dict):
            for kind in self.radii:
                self.atom_kinds[kind]['radius'] *= self.radii[kind]*scale
        #------------------------------------------------------------
        # bond_kinds
        if not self.bond_list and self.bond_cutoff:
            self.bond_list = get_bondpairs(self.atoms, cutoff = self.bond_cutoff)
        self.bond_kinds = get_bond_kinds(self.atoms, self.atom_kinds, self.bond_list)
        if self.debug:
            print(self.bond_list)
            print(self.bond_kinds)
        #------------------------------------------------------------
        self.polyhedra_kinds = get_polyhedra_kinds(self.atoms, self.atom_kinds, self.bond_list, polyhedra_dict = self.polyhedra_dict)
        #------------------------------------------------------------
        # cell 
        # disp = atoms.get_celldisp().flatten()
        if self.show_unit_cell == 'default':
            if self.atoms.pbc.any():
                self.show_unit_cell = True
            else:
                self.show_unit_cell = False
        self.cell_vertices = None
        if self.show_unit_cell:
            cell_vertices = np.empty((2, 2, 2, 3))
            for c1 in range(2):
                for c2 in range(2):
                    for c3 in range(2):
                        cell_vertices[c1, c2, c3] = np.dot([c1, c2, c3],
                                                           self.cell)
            cell_vertices.shape = (8, 3)
            self.cell_vertices = cell_vertices
        #
        if self.bbox is None:
            bbox = np.zeros([3, 2])
            # print(self.positions)
            R = self.atom_kinds[max(self.atom_kinds, key = lambda x: self.atom_kinds[x]['radius'])]['radius']
            for i in range(3):
                P1 = (self.positions[:, i] - R).min(0)
                P2 = (self.positions[:, i] + R).max(0)
                if self.show_unit_cell:
                    C1 = (self.cell_vertices[:, i] - self.celllinewidth).min(0)
                    C2 = (self.cell_vertices[:, i] + self.celllinewidth).max(0)
                    P1 = min(P1, C1)
                    P2 = max(P2, C2)
                bbox[i] = [P1, P2]
            self.bbox = bbox
            self.w = self.bbox[0][1] - self.bbox[0][0]
            self.h = self.bbox[1][1] - self.bbox[1][0]
        else:
            self.w = (self.bbox[0][1] - self.bbox[0][0]) * scale
            self.h = (self.bbox[1][1] - self.bbox[1][0]) * scale
        # print(self.bbox)
        # print(self.h, self.w)
        self.com = np.mean(self.bbox, axis=1)
        if self.camera_target is None:
            self.camera_target = self.com
        
        if not self.ortho_scale:
            if self.w > self.h:
                self.ortho_scale = self.w + 1
            else:
                self.ortho_scale = self.h + 1
        print(self.w, self.h, self.ortho_scale)
        #
        constr = self.atoms.constraints
        self.constrainatoms = []
        for c in constr:
            if isinstance(c, FixAtoms):
                for n, i in enumerate(c.index):
                    self.constrainatoms += [i]
        #
        self.material_styles_dict = {
            'jmol'    : {'Specular': 1.0, 'Roughness': 0.001, 'Metallic': 1.0},
            'ase3'    : {'Metallic': 1.0, 'Roughness': 0.001},
            'ceramic' : {'Subsurface': 0.1, 'Metallic': 0.02, 'Specular': 0.5, 'Roughness': 0.0},
            'plastic' : {'Metallic': 0.0, 'Specular': 0.5, 'Roughness': 0.7, 'Sheen Tint': 0.5, 'Clearcoat Roughness': 0.03, 'IOR': 1.6},
            'glass'   : {'Metallic': 0.0, 'Specular': 0.5, 'Roughness': 0.0, 'Clearcoat': 0.5, 'Clearcoat Roughness': 0.03, 'IOR': 1.45, 'Transmission': 0.98},
            'blase'   : {'Metallic': 0.02, 'Specular': 0.2, 'Roughness': 0.4, },
            'mirror'  : {'Metallic': 0.99, 'Specular': 2.0, 'Roughness': 0.001},
            }
        # ------------------------------------------------------------------------
        # remove all objects, Select objects by type
        clean_default()
        # 'AtomProp'.
        ALL_FRAMES = []
        # A list of ALL balls which are put into the scene
        self.STRUCTURE = []
        # COLLECTION
        # Before we start to draw the atoms, we first create a collection for the
        # atomic structure. All atoms (balls) are put into this collection.
        if self.build_collection:
            self.coll_name = os.path.basename(self.name) + '_blase'
            self.scene = bpy.context.scene
            self.coll = bpy.data.collections.new(self.coll_name)
            self.scene.collection.children.link(self.coll)
            self.coll_cell = bpy.data.collections.new('cell')
            self.coll_atom_kinds = bpy.data.collections.new('atoms')
            self.coll_bond_kinds = bpy.data.collections.new('bonds')
            self.coll_polyhedra_kinds = bpy.data.collections.new('polyhedras')
            self.coll_isosurface = bpy.data.collections.new('isosurfaces')
            self.coll_instancer = bpy.data.collections.new('instancers')
            self.coll.children.link(self.coll_cell)
            self.coll.children.link(self.coll_atom_kinds)
            self.coll.children.link(self.coll_bond_kinds)
            self.coll.children.link(self.coll_polyhedra_kinds)
            self.coll.children.link(self.coll_isosurface)
            self.coll.children.link(self.coll_instancer)
        # ------------------------------------------------------------------------
        # DRAWING THE ATOMS
        bpy.ops.object.select_all(action='DESELECT')
        # CAMERA and LIGHT SOURCES
        if self.camera:
            self.add_camera()
        if self.light:
            self.add_light()
        #
        if self.world:
            world = self.scene.world
            world.use_nodes = True
            node_tree = world.node_tree
            rgb_node = node_tree.nodes.new(type="ShaderNodeRGB")
            rgb_node.outputs["Color"].default_value = (1, 1, 1, 1)
            node_tree.nodes["Background"].inputs["Strength"].default_value = 1.0
            node_tree.links.new(rgb_node.outputs["Color"], node_tree.nodes["Background"].inputs["Color"])
        #
    def draw(self, coll = None):
        if not coll:
     	   coll = self.coll
        draw_cell(self, coll = coll)
        draw_atoms(self, coll = coll)
        # draw_bonds_2(self, coll = coll)
        draw_bonds(self, coll = coll)
        draw_polyhedras(self, coll = coll)
        if self.isosurface:
            volume = self.isosurface[0]
            icolor = 0
            for level in self.isosurface[1:]:
                draw_isosurface(self, coll = coll, volume=volume, level=level, icolor = icolor)
                icolor += 1

    #
    def look_at(self, obj, target, roll=0):
        """
        Rotate obj to look at target
        """
        if not isinstance(target, Vector):
            target = Vector(target)
        loc = obj.location
        direction = target - loc
        quat = direction.to_track_quat('-Z', 'Y')
        quat = quat.to_matrix().to_4x4()
        rollMatrix = Matrix.Rotation(roll, 4, 'Z')
        loc = loc.to_tuple()
        obj.matrix_world = quat @ rollMatrix
        obj.location = loc

    def add_camera(self):
        '''
        Type:   enum in ['PERSP', 'ORTHO'], default ORTHO
        '''
        # Create the camera
        target = self.camera_target
        camera_data = bpy.data.cameras.new("Camera")
        camera_data.lens = self.camera_lens
        camera_data.dof.aperture_fstop = self.fstop
        camera_data.type = self.camera_type
        camera = bpy.data.objects.new("Camera", camera_data)
        if self.camera_type == 'ORTHO' and not self.camera_loc:
            camera.location = [target[0], target[1], 200]
        else:
            camera.location = Vector(self.camera_loc)
        bpy.data.cameras['Camera'].type = self.camera_type
        if self.ortho_scale and self.camera_type == 'ORTHO':
            print('add_camera: ', self.ortho_scale)
            bpy.data.cameras['Camera'].ortho_scale = self.ortho_scale
        bpy.context.collection.objects.link(camera)
        self.look_at(camera, target, roll = radians(0))
        self.STRUCTURE.append(camera)
        bpy.context.scene.camera = camera
    def add_light(self, light_type = 'SUN', name = "Light", energy = 5):
        # Create the lamp
        '''
        POINT Point, Omnidirectional point light source.
        SUN Sun, Constant direction parallel ray light source.
        SPOT Spot, Directional cone light source.
        AREA Area, Directional area light source.
        '''
        if light_type:
            self.light_type = light_type
        light_data = bpy.data.lights.new(name=name, type=self.light_type)
        light_data.energy = energy
        lamp = bpy.data.objects.new(name, light_data)
        lamp.location = Vector(self.light_loc)
        bpy.context.collection.objects.link(lamp)
        # Some properties for cycles
        lamp.data.use_nodes = True
        lamp.data.node_tree.nodes['Emission'].inputs['Strength'].default_value = 0.1
        self.look_at(lamp, self.camera_target, roll = radians(0))
        self.STRUCTURE.append(lamp)

    
    def highlight_atoms(self, atomlist, shape = 'sphere', radius_scale=1.4,
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
        for index in atomlist:
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
            self.STRUCTURE.append(ball)
            coll_highlight.objects.link(ball)
    def draw_plane(self, color = (0.2, 0.2, 1.0, 1.0), size = 200, location = (0, 0, -1.0), bsdf_inputs = None, rotation = (0, 0, 0), material_style = 'blase'):
        """
        """
        # build materials
        if not bsdf_inputs:
            bsdf_inputs = self.material_styles_dict[material_style]
        material = bpy.data.materials.new('plane')
        material.name = 'plane'
        material.diffuse_color = color
        # material.blend_method = 'BLEND'
        material.use_nodes = True
        principled_node = material.node_tree.nodes['Principled BSDF']
        principled_node.inputs['Alpha'].default_value = color[3]
        for key, value in bsdf_inputs.items():
                principled_node.inputs[key].default_value = value
        # Instantiate a floor plane
        bpy.ops.mesh.primitive_plane_add(size=size, location=location, rotation=rotation)
        current_object = bpy.context.object
        current_object.data.materials.append(material)


    def render(self, outfile = None):
        
        # ------------------------------------------------------------------------
        # render settings
        if outfile:
            self.outfile = outfile
        self.directory = os.path.split(self.outfile)[0]
        if self.directory and not os.path.exists(self.directory):
                os.makedirs(self.directory)  # cp2k expects dirs to exist
        self.scene.render.image_settings.file_format = 'PNG'
        self.scene.render.engine = self.engine
        if self.engine.upper() == 'CYCLES' and self.gpu:
            self.scene.cycles.device = 'GPU'
            prefs = bpy.context.preferences.addons['cycles'].preferences
            # print(prefs.get_devices())
            for device in prefs.devices:
                # print(device)
                device.use = True
            self.scene.render.tile_x = 256
            self.scene.render.tile_y = 256
            self.scene.cycles.samples = self.num_samples

        self.scene.render.film_transparent = True
        # self.scene.render.alpha = 'SKY' # in ['TRANSPARENT', 'SKY']
        self.scene.render.resolution_x = self.resolution_x
        self.scene.render.resolution_y = int(self.resolution_x*self.h/self.w)
        print('dimension: ', self.scene.render.resolution_x, self.scene.render.resolution_y)
        self.scene.render.filepath = '{0}'.format(self.outfile)
        if self.save_to_blend:
            print('saving to {0}.blend'.format(self.outfile))
            bpy.ops.wm.save_as_mainfile('EXEC_SCREEN', filepath = '{0}.blend'.format(self.outfile))
        elif self.run_render:
            bpy.ops.render.render(write_still = 1, animation = self.animation)
    def export(self, filename = 'blender-ase.obj'):
        # render settings
        if filename.split('.')[-1] == 'obj':
            bpy.ops.export_scene.obj(filename)
        if filename.split('.')[-1] == 'x3d':
            bpy.ops.export_scene.x3d(filename)
        if filename.split('.')[-1] == 'xyz':
            self.export_xyz(filename)

    def render_move_camera(self, filename, loc1, loc2, n):
        # render settings
        from blase.tools import getEquidistantPoints
        locs = getEquidistantPoints(loc1, loc2, n)
        i = 0
        for loc in locs:
            bpy.data.objects['Camera'].location = loc
            self.render(self, filename + '_{0:3d}'.format(i))
            i += 1
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
            # for kind, datas in bond_kinds.items():
            #     for ibin in range(self.nbins):
            #         if len(datas['bins'][ibin]['vertices']) == 0:
            #             continue
            #         obj_bond = bpy.data.objects['bond_kind_{0}_{1}'.format(kind, ibin)]
            #         nverts = len(obj_bond.data.vertices)
            #         for j in range(nverts):
            #             obj_bond.data.vertices[j].co = datas['bins'][ibin]['vertices'][j]
            #             obj_bond.data.vertices[j].keyframe_insert('co', frame=i)
        self.scene.frame_start = 1
        self.scene.frame_end = self.nimages
    def export_xyz(self, filename):
        '''
        '''
        atoms = Atoms()
        for kind, datas in self.atom_kinds.items():
            obj_atom = bpy.data.objects['atom_kind_{0}'.format(kind)]
            nverts = len(obj_atom.data.vertices)
            for j in range(nverts):
                atom = Atom(self.atom_kinds[kind]['element'], obj_atom.data.vertices[j].co)
                atoms.append(atom)
        atoms.write(filename)

class ObjAse(bpy.types.PropertyGroup):
    element: bpy.props.StringProperty()
    index: bpy.props.IntProperty()
    x: bpy.props.FloatProperty()
    y: bpy.props.FloatProperty()
    z: bpy.props.FloatProperty()