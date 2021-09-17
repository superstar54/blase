"""Python module for drawing and rendering ase atoms objects using blender.
Basic idea:
"""
import bpy
from mathutils import Vector, Matrix
import os
import numpy as np
from math import pi, sqrt, radians, acos, atan2
from blase.tools import get_bbox
from blase.data import default_blase_settings, material_styles_dict
import logging
import sys

logging.basicConfig(stream=sys.stdout,
                    format=('%(levelname)-8s '
                            '[%(funcName)s]: %(message)s'),
                    level=logging.INFO)

logger = logging.getLogger('blase')


class Blase():
    """
    Blase object to render atomic structure.
    """
    #
    def __init__(self, output_image = 'bout', debug = False,
                              **parameters):
        for k, v in default_blase_settings.items():
            setattr(self, k, parameters.pop(k, v))
        #
        self.logger = logger
        if debug:
            print('debug is: ', debug)
            self.logger.setLevel(debug)
        #
        if 'Cube' in bpy.data.objects:
            bpy.data.objects.remove(bpy.data.objects["Cube"], do_unlink=True)
        self.output_image = output_image
        self.scene = bpy.context.scene
        if 'blase' in bpy.data.collections:
            self.coll = bpy.data.collections['blase']
        else:
            self.coll = bpy.data.collections.new('blase')
            self.scene.collection.children.link(self.coll)
        # ------------------------------------------------------------------------
        # CAMERA and LIGHT SOURCES
        if self.camera:
            self.logger.debug('Add camera')
            self.set_camera()
        if self.light:
            self.logger.debug('Add light')
            self.set_light()
        #
        if self.world:
            world = self.scene.world
            world.use_nodes = True
            node_tree = world.node_tree
            rgb_node = node_tree.nodes.new(type="ShaderNodeRGB")
            rgb_node.outputs["Color"].default_value = (1, 1, 1, 1)
            node_tree.nodes["Background"].inputs["Strength"].default_value = 1.0
            node_tree.links.new(rgb_node.outputs["Color"], node_tree.nodes["Background"].inputs["Color"])
    def set_camera(self, camera_type = None, camera_lens = None,):
        '''
        Set camera.

        camera_type: str

            * PERSP
            * ORTHO

        camera_lens: float
        '''
        # check camera exist or not
        if not camera_type:
            camera_type = self.camera_type
        if not camera_lens:
            camera_lens = self.camera_lens
        if 'Camera_blase' in bpy.data.objects:
            camera = bpy.data.objects['Camera_blase']
        else:
            camera_data = bpy.data.cameras.new("Camera_blase")
            camera = bpy.data.objects.new("Camera_blase", camera_data)
            self.coll.objects.link(camera)
        camera.data.lens = camera_lens
        camera.data.dof.aperture_fstop = self.fstop
        camera.data.type = camera_type
        # image size and target
        # print('bbox: ', self.bbox)
        self.width = (self.bbox[0][1] - self.bbox[0][0])
        self.height = (self.bbox[1][1] - self.bbox[1][0])
        self.com = np.mean(self.bbox, axis=1)
        if self.camera_target is None:
            self.camera_target = self.com
        if not self.ortho_scale:
            if self.width > self.height:
                self.ortho_scale = self.width + 1
            else:
                self.ortho_scale = self.height + 1
        print(self.width, self.height)
        print('ortho_scale: ', self.ortho_scale)
        target = self.camera_target
        if camera_type == 'ORTHO' and not self.camera_loc:
            camera.location = [target[0], target[1], 200]
        else:
            camera.location = Vector(self.camera_loc)
        if self.ortho_scale and camera_type == 'ORTHO':
            print('set_camera: ', self.ortho_scale)
            camera.data.ortho_scale = self.ortho_scale
        self.look_at(camera, target, roll = radians(0))
        bpy.context.scene.camera = camera
    def set_light(self, location = None, light_type = 'SUN', name = "Light", energy = 10):
        '''
        location: array

        light_type: str

            * POINT, Omnidirectional point light source.
            * SUN, Constant direction parallel ray light source.
            * SPOT, Directional cone light source.
            * AREA, Directional area light source.
        
        energy: float

        '''
        # check light exist or not
        if 'Light_blase' in bpy.data.objects:
            light = bpy.data.objects['Light_blase']
        else:
            light_data = bpy.data.lights.new("Light_blase", type = light_type)
            light = bpy.data.objects.new("Light_blase", light_data)
            self.coll.objects.link(light)
        light.data.type = light_type
        light.data.energy = energy
        if location:
            light.location = Vector(location)
        else:
            light.location = Vector(self.light_loc)
        # Some properties for cycles
        light.data.use_nodes = True
        light.data.node_tree.nodes['Emission'].inputs['Strength'].default_value = 0.1
        self.look_at(light, self.camera_target, roll = radians(0))
        print('set light: ', self.ortho_scale)
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
    def draw_plane(self, location = (0, 0, -1.0), color = (0.2, 0.2, 1.0, 1.0), size = 200, bsdf_inputs = None, material_style = 'blase'):
        """
        Draw a plane.

        location: array

        color: array

        size: float

        bsdf_inputs: dict

        material_style: str

        """
        # build materials
        if not bsdf_inputs:
            bsdf_inputs = material_styles_dict[material_style]
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

    def render(self, output_image = None):
        """
        """
        # render settings
        if output_image:
            self.output_image = output_image
        self.directory = os.path.split(self.output_image)[0]
        if self.directory and not os.path.exists(self.directory):
                os.makedirs(self.directory)  # cp2k expects dirs to exist
        self.scene.render.image_settings.file_format = 'PNG'
        self.scene.render.engine = self.engine
        if self.engine.upper() == 'BLENDER_WORKBENCH' and 'StudioLight_blase.sl' in bpy.context.preferences.studio_lights:
            bpy.data.scenes['Scene'].display.shading.studio_light = 'StudioLight_blase.sl'
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
        self.scene.render.resolution_y = int(self.resolution_x*self.height/self.width)
        print('dimension: ', self.scene.render.resolution_x, self.scene.render.resolution_y)
        self.scene.render.filepath = '{0}'.format(self.output_image)
        if self.save_to_blend:
            print('saving to {0}.blend'.format(self.output_image))
            bpy.ops.wm.save_as_mainfile('EXEC_SCREEN', filepath = '{0}.blend'.format(self.output_image))
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
