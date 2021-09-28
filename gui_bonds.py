import bpy
import bmesh
from mathutils import Vector
from bpy.types import (Panel,
                       Operator,
                       AddonPreferences,
                       PropertyGroup,
                       )
from bpy.props import (StringProperty,
                       BoolProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       EnumProperty,
                       PointerProperty,
                       )

from blase.gui_io import import_blase


import os
from math import sqrt
from copy import copy
from ase import Atom, Atoms
from ase.build import molecule, bulk
import json
from blase.render import Blase
from blase.butils import read_batoms_collection_list, read_atoms_list


# The panel.
class Bonds_PT_prepare(Panel):
    bl_label       = "Blase Bonds"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    # bl_options     = {'DEFAULT_CLOSED'}
    bl_category = "Bonds"
    bl_idname = "BLASE_PT_Bonds"

  
    def draw(self, context):
        layout = self.layout
        bopanel = context.scene.bopanel

        box = layout.box()
        row = box.row()
        row.prop(bopanel, "collection_list")
        row = box.row()
        row.prop(bopanel, "atoms_list")

        box = layout.box()
        # row = box.row()
        # row.label(text="Materials")
        row = box.row()
        row.prop(bopanel, "materials")


        row = box.row()
        row.prop(bopanel, "radius")
       
        col = box.column(align=True)
        col.label(text="Color")
        row = box.row()
        row.prop(bopanel, "colorR")
        row.prop(bopanel, "colorG")
        row.prop(bopanel, "colorB")



class BondsProperties(bpy.types.PropertyGroup):
    def Callback_collection_list(self, context):
        items = read_batoms_collection_list()
        items = [(item, item, "") for item in items]
        items = tuple(items)
        return items
    def Callback_materials(self, context):
        bopanel = bpy.context.scene.bopanel
        modify_materials(bopanel.atoms_list, bopanel.materials)
    def Callback_atoms_list(self, context):
        bopanel = bpy.context.scene.bopanel
        coll = bpy.data.collections[bopanel.collection_list]
        items = read_atoms_list(coll)
        items = [(item, item, "") for item in items]
        items = tuple(items)
        return items
    collection_list: EnumProperty(
        name="Collection",
        description="Collection",
        items=Callback_collection_list)
    atoms_list: EnumProperty(
        name="Atom",
        description="Atom",
        items=Callback_atoms_list)
    materials: EnumProperty(
        name="Materials",
        description="Materials models",
        items=(('0',"blase", ""),
               ('1',"metal", ""),
               ('2',"glass",""),
               ('3',"plastic", ""),
               ('4',"mirror", "")),
               default='0', update=Callback_materials)

    radius: FloatProperty(
        name="Radius", default=False,
        description = "Radius")
    colorR: FloatProperty(
        name="R", default=False,
        description = "R")
    colorG: FloatProperty(
        name="G", default=False,
        description = "G")
    colorB: FloatProperty(
        name="B", default=False,
        description = "B")
    



# Modifying the radius of a selected atom or stick
def modify_materials(collection_name, model_type, atoms = None):
    # Modify atom radius (all selected)
    pass