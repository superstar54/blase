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
from blase.bio import Blase
from blase.btools import draw_cell, draw_atoms, draw_bonds, draw_polyhedras, draw_isosurface, bond_source, cylinder_mesh_from_instance, clean_default
from blase.utils import read_blase_collection_list, read_blase_collection, read_atoms_list


# The panel.
class Atoms_PT_prepare(Panel):
    bl_label       = "Blase Atoms"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    # bl_options     = {'DEFAULT_CLOSED'}
    bl_category = "Atoms"
    bl_idname = "BLASE_PT_Atoms"

  
    def draw(self, context):
        layout = self.layout
        atpanel = context.scene.atpanel

        box = layout.box()
        row = box.row()
        row.prop(atpanel, "collection_list")
        row = box.row()
        row.prop(atpanel, "atoms_list")

        row = box.row()
        row.prop(atpanel, "materials")
        row = box.row()
        row.prop(atpanel, "radius")
       
        col = box.column(align=True)
        row = box.row()
        row.prop(atpanel, "color")



class AtomsProperties(bpy.types.PropertyGroup):
    def Callback_collection_list(self, context):
        items = read_blase_collection_list()
        items = [(item, item, "") for item in items]
        items = tuple(items)
        return items
    def Callback_materials(self, context):
        atpanel = bpy.context.scene.atpanel
        modify_materials(atpanel.atoms_list, atpanel.materials)
    def Callback_atoms_list(self, context):
        atpanel = bpy.context.scene.atpanel
        coll = bpy.data.collections[atpanel.collection_list]
        items = read_atoms_list(coll)
        items = [(item, item, "") for item in items]
        items = tuple(items)
        return items
    def Callback_modify_radius(self, context):
        atpanel = bpy.context.scene.atpanel
        print('Callback_modify_radius')
        modify_radius(atpanel.collection_list, atpanel.atoms_list, atpanel.radius)
    def Callback_modify_color(self, context):
        atpanel = bpy.context.scene.atpanel
        print('Callback_modify_radius')
        color = atpanel.color
        modify_color(atpanel.collection_list, atpanel.atoms_list, color)


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
        name="Radius", default=1.0,
        description = "Radius", update = Callback_modify_radius)
    colorR: FloatProperty(
        name="R", default=False,
        description = "R")
    colorG: FloatProperty(
        name="G", default=False,
        description = "G")
    colorB: FloatProperty(
        name="B", default=False,
        description = "B", update = Callback_modify_color)
    color = FloatVectorProperty(  
        name="color",
        subtype='COLOR',
        default=(1.0, 1.0, 1.0),
        min=0.0, max=1.0,
        description="color picker"
        )
    



# Modifying the radius of a selected atom or stick
def modify_materials(collection_name, model_type, atoms = None):
    coll = bpy.data.collections[collection_name]
    if not batoms:
        batoms = read_blase_collection(coll)
    print('drawing atoms')
    batoms.draw_atoms()

# Modifying the radius of a selected atom or stick
def modify_radius(collection_name, atoms_list, radius, batoms = None):
    # Modify atom radius (all selected)
    coll = bpy.data.collections[collection_name]
    if not batoms:
        batoms = read_blase_collection(coll)
    print('drawing atoms')
    batoms.draw_atoms(props={atoms_list:{'radius': radius,}})
# Modifying the radius of a selected atom or stick
def modify_color(collection_name, atoms_list, color, batoms = None):
    # Modify atom radius (all selected)
    coll = bpy.data.collections[collection_name]
    if not batoms:
        batoms = read_blase_collection(coll)
    print('drawing atoms')
    batoms.draw_atoms(props={atoms_list:{'color': color,}})