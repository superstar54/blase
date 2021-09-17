import json
import bpy
from bpy.types import Operator, AddonPreferences
from bpy_extras.io_utils import ImportHelper, ExportHelper
from bpy.props import (
        StringProperty,
        BoolProperty,
        EnumProperty,
        IntProperty,
        FloatProperty,
        )
from ase.io.cube import read_cube_data
import pickle
from blase.bio import Blase
from blase.batoms import Batoms

# -----------------------------------------------------------------------------
#                                                                     Operators

# This is the class for the file dialog.
class IMPORT_OT_blase(Operator, ImportHelper):
    bl_idname = "import_mesh.blase"
    bl_label  = "Import blase (*.blase)"
    bl_options = {"PRESET", "REGISTER", "UNDO"}

    filename_ext = ".blase"
    # filter_glob: StringProperty(default="*.blase", options={'HIDDEN'},)

    camera: BoolProperty(
        name="Camera", default=False,
        description="Do you need a camera?")
    light: BoolProperty(
        name="Light", default=False,
        description = "Do you need a light?")
    world: BoolProperty(
        name="World", default=False,
        description = "Do you need a world light?")
    model_type: EnumProperty(
        name="Type",
        description="Choose model",
        items=(('0',"Space-filling", "Use ball"),
               ('1',"Ball-and-stick", "Use ball and stick"),
               ('2',"Polyhedral","Use polyhedral"),
               ('3',"Stick", "Use stick")),
               default='0',)
    ball_type: EnumProperty(
        name="Type of ball",
        description="Choose model",
        items=(('0',"Ball", "Use ball and stick"),
               ('1',"Meta", "Use ball")),
               default='0',)
    label: StringProperty(
        name = "label", 
        description = "Label")
    bond_cutoff: FloatProperty(
        name = "Bond_cutoff", default=1.0, min=0.0001,
        description = "Bond cutoff")
    show_unit_cell: BoolProperty(
        name = "Show_unit_cell", default=False,
        description = "Show unit cell")
    atomradius: EnumProperty(
        name="Type of radius",
        description="Choose type of atom radius",
        items=(('0', "Pre-defined", "Use pre-defined radius"),
               ('1', "Atomic", "Use atomic radius"),
               ('2', "van der Waals", "Use van der Waals radius")),
               default='0',)
    make_real: FloatProperty(
        name = "Make_real", default=False,
        description = "")


    def draw(self, context):
        
        layout = self.layout
        box = layout.box()
        row = box.row()
        row.label(text="Adding Structure")
        box = layout.box()
        row = box.row()
        row.prop(self, "label")
        #
        row = layout.row()
        row.prop(self, "camera")
        row.prop(self, "light")
        row.prop(self, "world")
        #
        box = layout.box()
        row = box.row()
        # row.label(text="Unit_cell")
        row.prop(self, "show_unit_cell")
        # Balls
        box = layout.box()
        row = box.row()
        row.label(text="Structure model")
        row = box.row()
        col = row.column()
        col.prop(self, "model_type")
        box = layout.box()
        row = box.row()
        row.active = (self.model_type == "1")
        #
        

    def execute(self, context):

        # This is to determine the path.
        self.inputfile = bpy.path.abspath(self.filepath)

        # Execute main routine
        import_blase(self.inputfile, 
               self.model_type,
               self.camera,
               self.light,
               self.world,
               self.show_unit_cell,
               self.label,
               )

        return {'FINISHED'}


def import_blase(inputfile, 
               model_type = '0',
               camera = 'True',
               light = 'True',
               world = 'False',
               show_unit_cell = False,
               label = None,
               ):
    #
    from ase.io import read
    print('='*30)
    print('Import structure')
    print('='*30)
    print(model_type)
    #
    kwargs = {'show_unit_cell': show_unit_cell, 
          'world': world,
          'camera': camera,
          'light': light,
          }
    if inputfile.split('.')[-1] == 'cube':
        images, data = read_cube_data(inputfile)
    else:
        images = read(inputfile)
    if not label:
        label = inputfile
    bobj = Batoms(label = label, atoms = images, model_type=model_type) #, **kwargs)
    bobj.draw()
