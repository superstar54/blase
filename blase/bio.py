# ###################################################################
"""Python module for drawing and rendering ase atoms objects using blender.

Part of code based on ase.io.pov, 
https://wiki.fysik.dtu.dk/ase/dev/_modules/ase/io/pov.html

Part of the utils from
 https://github.com/yuki-koyama/blender-cli-rendering

"""

# ###################################################################

import bpy
from mathutils import Vector, Matrix
import os
import numpy as np
from ase.io.utils import PlottingVariables, cell_to_lines
from ase.constraints import FixAtoms
from ase.utils import basestring
from ase.io import read, write
from ase.build import molecule, bulk
from ase.data import covalent_radii, atomic_numbers
from ase.data.colors import jmol_colors
from math import pi, sqrt, radians, acos, atan2
from skimage import measure
from blase.tools import get_atom_kinds, get_bond_kinds, get_bondpairs
from blase.btools import bond_source, bond_mesh_from_instance, removeAll


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
        'camera_loc': [0, 0, 50],  # x, y is the image plane, z is *out* of the screen
        'camera_type': 'ORTHO',  #  ['PERSP', 'ORTHO']
        'ortho_scale': None, #
        'camera_lens': 10,  #
        'fstop': 0.5,
        'camera_target': None, #
        'world': False,
        'light': True,
        'light_type': 'SUN', # 'POINT', 'SUN', 'SPOT', 'AREA'
        'point_lights': [],  # 
        'light_strength': 1.0,
        'background': 'White',  # color
        'textures': None,  # length of atoms list of texture names
        'engine': 'BLENDER_EEVEE', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH'
        'transmits': None,  # transmittance of the atoms
        'show_unit_cell': 0,
        'celllinewidth': 0.05,  # radius of the cylinders representing the cell
        'bbox': None,
        'bondlinewidth': 0.10,  # radius of the cylinders representing bonds
        'balltypes': None,
        'radii': None, 
        'colors': None,
        'kind_props': None,
        'bond_cutoff': None,  # 
        'bondlist': [],  # [[atom1, atom2], ... ] pairs of bonding atoms
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
        'outfile': 'bout',
        }  

    def __init__(self, images, rotations=None, scale=1,
                              **parameters):
        for k, v in self.default_settings.items():
            setattr(self, k, parameters.pop(k, v))
        #
        if not isinstance(images, list):
            images = [images]
        self.nimages = len(images)
        if rotations:
            for rotation in rotations:
                for i in range(self.nimages):
                    images[i].rotate(rotation[0], rotation[1], rotate_cell = True)
        #
        self.images = images
        self.atoms = images[0]
        self.natoms = len(self.atoms)
        self.positions = self.atoms.positions*scale
        self.numbers = self.atoms.get_atomic_numbers()
        self.symbols = self.atoms.get_chemical_symbols()
        self.cell = self.atoms.get_cell()*scale
        self.celllinewidth = self.celllinewidth*scale
        self.bondlinewidth = self.bondlinewidth*scale
        self.atom_kinds = get_atom_kinds(self.atoms, self.kind_props)
        self.nkinds = len(self.atom_kinds)
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
        
        # disp = atoms.get_celldisp().flatten()
        self.cell_vertices = None
        if self.show_unit_cell > 0:
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
                if self.show_unit_cell > 0:
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
                # self.resolution_y = 
            else:
                self.ortho_scale = self.h + 1
        # print(self.bbox)
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
        removeAll()
        # 'AtomProp'.
        ALL_FRAMES = []
        # A list of ALL balls which are put into the scene
        self.STRUCTURE = []
        # COLLECTION
        # Before we start to draw the atoms, we first create a collection for the
        # atomic structure. All atoms (balls) are put into this collection.
        filename = 'ase'
        self.coll_structure_name = os.path.basename(filename)
        self.scene = bpy.context.scene
        self.coll_structure = bpy.data.collections.new(self.coll_structure_name)
        self.scene.collection.children.link(self.coll_structure)
        self.coll_cell = bpy.data.collections.new('cell')
        self.coll_atom_kinds = bpy.data.collections.new('atoms')
        self.coll_bond_kinds = bpy.data.collections.new('bonds')
        self.coll_isosurface = bpy.data.collections.new('isosurfaces')
        self.coll_structure.children.link(self.coll_cell)
        self.coll_structure.children.link(self.coll_atom_kinds)
        self.coll_structure.children.link(self.coll_bond_kinds)
        self.coll_structure.children.link(self.coll_isosurface)
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
        camera = bpy.data.objects.new("Camera", camera_data)
        camera.location = Vector(self.camera_loc)
        if self.camera_type == 'ORTHO':
            camera.location = [target[0], target[1], self.camera_loc[2]]
        bpy.data.cameras['Camera'].type = self.camera_type
        if self.ortho_scale and self.camera_type == 'ORTHO':
            bpy.data.cameras['Camera'].ortho_scale = self.ortho_scale
        bpy.context.collection.objects.link(camera)
        self.look_at(camera, target, roll = radians(0))
        self.STRUCTURE.append(camera)
        bpy.context.scene.camera = camera
    def add_light(self, loc = [0, 0, 200], light_type = 'SUN', name = "Light", energy = 5):
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
        lamp.location = Vector(loc)
        bpy.context.collection.objects.link(lamp)
        # Some properties for cycles
        lamp.data.use_nodes = True
        lamp.data.node_tree.nodes['Emission'].inputs['Strength'].default_value = 0.1
        self.STRUCTURE.append(lamp)
    #
    def draw_cylinder(self, p1, p2, material, radius = 0.1):
        """
        Draw cylinder connect two points.
        """
        diff = tuple([c2-c1 for c2, c1 in zip(p1, p2)])
        center = tuple([(c2+c1)/2 for c2, c1 in zip(p1, p2)])
        magnitude = pow(sum([(c2-c1)**2
                        for c1, c2 in zip(p1, p2)]), 0.5)
        # Euler rotation calculation, (Vector from mathutils, acos from math)
        bpy.ops.mesh.primitive_cylinder_add(radius=radius, depth=magnitude, location=center)
        phi = atan2(diff[1], diff[0])
        theta = acos(diff[2]/magnitude)
        bpy.context.object.rotation_euler[1] = theta
        bpy.context.object.rotation_euler[2] = phi
        cylinder = bpy.context.view_layer.objects.active
        cylinder.active_material = material
        return cylinder
    #========================================================
    def draw_cell(self, ):
        """
        Draw unit cell
        """
        if self.cell_vertices is not None:
            self.cell_vertices.shape = (2, 2, 2, 3)
            # build materials
            material = bpy.data.materials.new('cell')
            material.name = 'cell'
            material.diffuse_color = (0.8, 0.25, 0.25, 1.0)
            #
            ic = 0
            for c in range(3):
                for j in ([0, 0], [1, 0], [1, 1], [0, 1]):
                    # print('cell: ', ic)
                    parts = []
                    for i in range(2):
                        j.insert(c, i)
                        parts.append(self.cell_vertices[tuple(j)])
                        del j[c]

                    distance = np.linalg.norm(parts[1] - parts[0])
                    if distance < 1e-12:
                        continue
                    #
                    cylinder = self.draw_cylinder(parts[0], parts[1], material, radius = self.celllinewidth)
                    self.STRUCTURE.append(cylinder)
                    self.coll_cell.objects.link(cylinder)
                    ic += 1
    #
    def draw_atoms(self, bsdf_inputs = None, material_style = 'blase'):
        '''
        Draw atoms
        bsdf_inputs: dict
            The key and value for principled_bsdf node
        material_style: string
            Select materials type from ['blase', 'glass', 'ceramic', 'plastic'].
        '''
        # build materials
        self.materials = []
        if not bsdf_inputs:
            bsdf_inputs = self.material_styles_dict[material_style]
        for kind, datas in self.atom_kinds.items():
            print('atoms: ', kind)
            material = bpy.data.materials.new('atom_kind_{0}'.format(kind))
            material.diffuse_color = np.append(datas['color'], datas['transmit'])
            # material.blend_method = 'BLEND'
            material.use_nodes = True
            principled_node = material.node_tree.nodes['Principled BSDF']
            principled_node.inputs['Alpha'].default_value = datas['transmit']
            for key, value in bsdf_inputs.items():
                principled_node.inputs[key].default_value = value
            datas['materials'] = material
            #
            # print(datas['positions'][0])
            if datas['balltype'] == 'meta':
                print('metaball', kind)
                mbdata = bpy.data.metaballs.new('atom_kind'.format(kind))
                mbdata.render_resolution = 0.1
                mbdata.resolution = 0.2
                obj_atom = bpy.data.objects.new('meta_atom_kind'.format(kind), mbdata)
                obj_atom.data.materials.append(material)
                for co in datas['positions']:
                    mbele = mbdata.elements.new(type = 'BALL')
                    mbele.co = co
                    mbele.radius = datas['radius']*2.0
                    mbele.stiffness = 1.0
                self.STRUCTURE.append(obj_atom)
                self.coll_atom_kinds.objects.link(obj_atom)
            else:
                bpy.ops.mesh.primitive_uv_sphere_add(radius = datas['radius']) #, segments=32, ring_count=16)
                ball = bpy.context.view_layer.objects.active
                ball.data.materials.append(material)
                mesh = bpy.data.meshes.new('mesh_kind_{0}'.format(kind) )
                obj_atom = bpy.data.objects.new('atom_kind_{0}'.format(kind), mesh )
                # Associate the vertices
                obj_atom.data.from_pydata(datas['positions'], [], [])
                bpy.context.view_layer.objects.active = obj_atom
                # Make the object parent of the cube
                bpy.ops.object.shade_smooth()
                ball.parent = obj_atom
                # Make the object dupliverts
                obj_atom.instance_type = 'VERTS'
                ball.hide_set(True)
                self.STRUCTURE.append(obj_atom)
                self.coll_atom_kinds.objects.link(obj_atom)

    def draw_bonds(self, bond_rcutoff= 1.0, bsdf_inputs = None, material_style = 'blase'):
        '''
        Draw atom bonds
        '''
       
        ib = 0
        if not bsdf_inputs:
            bsdf_inputs = self.material_styles_dict[material_style]
        if self.bond_cutoff:
            self.bondlist = get_bondpairs(self.atoms, cutoff = self.bond_cutoff)
        #
        bond_kinds = get_bond_kinds(self.atoms, self.bondlist)
        # import pprint
        # pprint.pprint(bond_kinds)
        source = bond_source()
        for kind, datas in bond_kinds.items():
            print('bonds: ', kind)
            material = bpy.data.materials.new('bond_kind_{0}'.format(kind))
            material.diffuse_color = np.append(self.atom_kinds[kind]['color'], self.atom_kinds[kind]['transmit'])
            # material.blend_method = 'BLEND'
            material.use_nodes = True
            principled_node = material.node_tree.nodes['Principled BSDF']
            principled_node.inputs['Alpha'].default_value = self.atom_kinds[kind]['transmit']
            for key, value in bsdf_inputs.items():
                principled_node.inputs[key].default_value = value
            datas['materials'] = material
            #
            # print(datas['normals'])
            verts, faces = bond_mesh_from_instance(datas['centers'], datas['normals'], datas['lengths'], self.bondlinewidth, source)
            # print(datas['vertices'][0])
            # bpy.ops.mesh.primitive_uv_sphere_add(radius = datas['radius'])
            # create new mesh structure
            mesh = bpy.data.meshes.new("mesh_kind_{0}".format(kind))
            mesh.from_pydata(verts, [], faces)  
            mesh.update()
            for f in mesh.polygons:
                f.use_smooth = True
            obj_bond = bpy.data.objects.new("bond_kind_{0}".format(kind), mesh)
            obj_bond.data = mesh
            obj_bond.data.materials.append(material)
            bpy.ops.object.shade_smooth()

            # mesh = bpy.data.meshes.new('mesh_kind_{0}'.format(kind) )
            # obj_bond = bpy.data.objects.new('bond_kind_{0}'.format(kind), mesh )
            # # Associate the vertices
            # obj_bond.data.from_pydata(datas['vertices'], [], datas['faces'])
            # # Link the object to the scene and set it active
            # bpy.context.view_layer.objects.active = obj_bond
            # Make the object parent of the cube
            # Make the object dupliverts
            self.STRUCTURE.append(obj_bond)
            self.coll_bond_kinds.objects.link(obj_bond)

    def draw_isosurface(self, volume, level,
                        closed_edges = False, gradient_direction = 'descent',
                        color=(0.85, 0.80, 0.25) , transmit=1.0,
                        verbose = False, step_size = 1, 
                        bsdf_inputs = None, material_style = 'blase'):
        """Computes an isosurface from a volume grid.
        
        Parameters: 

        """
        cell = self.cell
        self.cell_vertices.shape = (2, 2, 2, 3)
        cell_origin = self.cell_vertices[0,0,0]
        #
        spacing = tuple(1.0/np.array(volume.shape))
        scaled_verts, faces, normals, values = measure.marching_cubes_lewiner(volume, level = level,
                        spacing=spacing,gradient_direction=gradient_direction , 
                        allow_degenerate = False, step_size=step_size)
        #
        scaled_verts = list(scaled_verts)
        nverts = len(scaled_verts)
        # transform
        for i in range(nverts):
            scaled_verts[i] = scaled_verts[i].dot(cell)
            scaled_verts[i] -= cell_origin
        faces = list(faces)
        print('Draw isosurface...')
        # print('verts: ', scaled_verts[0:5])
        # print('faces: ', faces[0:5])
        #material
        if not bsdf_inputs:
            bsdf_inputs = self.material_styles_dict[material_style]
        material = bpy.data.materials.new('isosurface')
        material.name = 'isosurface'
        material.diffuse_color = color + (transmit,)
        # material.alpha_threshold = 0.2
        # material.blend_method = 'BLEND'
        material.use_nodes = True
        principled_node = material.node_tree.nodes['Principled BSDF']
        principled_node.inputs['Alpha'].default_value = transmit
        for key, value in bsdf_inputs.items():
                principled_node.inputs[key].default_value = value
        #
        # create new mesh structure
        isosurface = bpy.data.meshes.new("isosurface")
        isosurface.from_pydata(scaled_verts, [], faces)  
        isosurface.update()
        for f in isosurface.polygons:
            f.use_smooth = True
        iso_object = bpy.data.objects.new("isosurface", isosurface)
        iso_object.data = isosurface
        iso_object.data.materials.append(material)
        bpy.ops.object.shade_smooth()
        self.coll_isosurface.objects.link(iso_object)
        self.STRUCTURE.append(iso_object)


    def highlight_atoms(self, atomlist, shape = 'sphere', radius_scale=1.4,
                           color=(0.0, 1.0, 0.0), transmit=0.6):
        """
        """
        # Draw atoms
        #
        coll_highlight = bpy.data.collections.new('highlight')
        self.coll_structure.children.link(coll_highlight)
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
            # bond_kinds = get_bond_kinds(self.images[i], self.bondlist, self.nbins)
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




