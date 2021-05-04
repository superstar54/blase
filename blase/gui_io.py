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
from ase import Atom, Atoms
from ase.io import read
from ase.io.cube import read_cube_data
import pickle
from blase.bio import Blase
from blase.btools import draw_cell, draw_atoms, draw_bonds, draw_polyhedras, draw_isosurface, bond_source, cylinder_mesh_from_instance, clean_default
from blase.batoms import Batoms

# -----------------------------------------------------------------------------
#                                                                     Operators

# This is the class for the file dialog.
class IMPORT_OT_blase(Operator, ImportHelper):
    bl_idname = "import_mesh.blase"
    bl_label  = "Import blase (*.blase)"
    bl_options = {'PRESET', 'UNDO'}

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
    scale: FloatProperty(
        name = "scale", default=1.0, min=0.0001,
        description = "Scale factor for all atom scale")
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
        box = layout.box()
        row = box.row()
        row.label(text="Scaling scale")
        box = layout.box()
        row = box.row()
        row.prop(self, "scale")
        box = layout.box()
        row = box.row()
        # row.label(text="Bond cutoff")
        row.prop(self, "bond_cutoff")
        
        box = layout.box()
        row = box.row()
        # row.label(text="Search_pbc")
        row.prop(self, "search_pbc")
        # Frames


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
               self.ball_type,
               self.scale,
               self.bond_cutoff,
               self.make_real,
               )

        return {'FINISHED'}


def import_blase(inputfile, 
               model_type = '0',
               camera = 'True',
               light = 'True',
               world = 'False',
               show_unit_cell = False,
               ball_type = 'Ball-and-stick',
               scale = 1.0,
               bond_cutoff = 1.0,
               search_pbc_atoms = False,
               search_molecule = False,
               make_real = False,
               name = None,
               ):
    #
    print('='*30)
    print('Import structure')
    print('='*30)
    print(model_type)
    #
    kwargs = {'show_unit_cell': show_unit_cell, 
          'scale': scale,
          'bond_cutoff': bond_cutoff,
          'world': world,
          'camera': camera,
          'light': light,
          'search_pbc_atoms':search_pbc_atoms,
          'search_molecule': search_molecule,
          'make_real': make_real,
          }
    if isinstance(inputfile, str):
        if inputfile.split('.')[-1] == 'cube':
            images, data = read_cube_data(inputfile)
        else:
            images = read(inputfile)
    else:
        images = inputfile
    print(images)
    bobj = Batoms(images, name = name, model_type=model_type) #, **kwargs)
    bobj.draw()








def export_blase():
    '''

    '''
    atoms = Atoms()
    cell_vertexs = []
    for obj in bpy.context.selected_objects:
        if "BOND" in obj.name.upper():
            continue
        if obj.type not in {'MESH', 'SURFACE', 'META'}:
            continue
        name = ""
        if 'atom_kind' in obj.name:
            print(obj.name)
            ele = obj.name[10:].split('.')[0]
            if len(obj.children) != 0:
                for vertex in obj.data.vertices:
                    location = obj.matrix_world @ vertex.co
                    atoms.append(Atom(ele, location))
            else:
                if not obj.parent:
                    location = obj.location
                    atoms.append(Atom(ele, location))
        # cell
        if 'point_cell' in obj.name:
            print(obj.name)
            if len(obj.children) != 0:
                for vertex in obj.data.vertices:
                    location = obj.matrix_world @ vertex.co
                    # print(location)
                    cell_vertexs.append(location)
    # print(atoms)
    # print(cell_vertexs)
    if cell_vertexs:
        cell = [cell_vertexs[4], cell_vertexs[2], cell_vertexs[1]]
        atoms.cell = cell
    return atoms

