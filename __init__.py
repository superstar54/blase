# Blase TOOLBAR  - Addon - Blender 2.9x
#
# THIS SCRIPT IS LICENSED UNDER GPL,
# please read the license block.

# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####


###
bl_info = {
    "name": "Blase toolbar",
    "author": "Xing Wang",
    "version": (0, 5),
    "blender": (2, 90, 0),
    "location": "File -> Import -> Blase (xyz, cif, pdb, ...)",
    "description": "Python module for drawing and rendering ASE (Atomic Simulation Environment) atoms and molecules objects using blender.",
    "warning": "",
    "category": "Import-Export",
}

from blase.batoms import Batoms
from blase.batom import Batom

import bpy
from bpy.types import (Panel,
                       Operator,
                       AddonPreferences,
                       PropertyGroup,
                       )
from . import (
        gui_io,  # Blase import/export
        gui_blase,  # Panel to edit atomic structure interactively.
        gui_bonds,  #
        gui_atoms,  #
        gui_cell,  #
        gui_uilist,
        )

def menu_func_import_blase(self, context):
    lay = self.layout
    lay.operator(gui_io.IMPORT_OT_blase.bl_idname,text="blase file (xyz, cif, pdb, ...)")


# Register
classes = [
        gui_io.IMPORT_OT_blase,
        gui_blase.Blase_PT_prepare,
        gui_blase.BlaseProperties,
        gui_blase.BlaseSettings,
        # gui_blase.BatomSettings,
        gui_blase.BlaseAtom,
        gui_blase.BlaseBond,
        gui_blase.ExportAtom,
        gui_blase.SplitAtom,
        gui_blase.AddMolecule,
        gui_blase.AddBulk,
        gui_blase.AddAtoms,
        gui_blase.CopyAtoms,
        gui_bonds.Bonds_PT_prepare,
        gui_bonds.BondsProperties,
        gui_atoms.Atoms_PT_prepare,
        gui_atoms.AtomsProperties,
        gui_cell.Cell_PT_prepare,
        gui_cell.CellProperties,
    ]
def register():
    bpy.types.TOPBAR_MT_file_import.append(menu_func_import_blase)
    for cls in classes:
        bpy.utils.register_class(cls)
    scene = bpy.types.Scene
    scene.blpanel = bpy.props.PointerProperty(type=gui_blase.BlaseProperties)
    scene.atpanel = bpy.props.PointerProperty(type=gui_atoms.AtomsProperties)
    scene.bopanel = bpy.props.PointerProperty(type=gui_atoms.AtomsProperties)
    scene.clpanel = bpy.props.PointerProperty(type=gui_cell.CellProperties)
    bpy.types.Collection.is_batoms = bpy.props.BoolProperty(name = 'is_batoms')
    bpy.types.Collection.blase = bpy.props.PointerProperty(name = 'blase', type = gui_blase.BlaseSettings)
    bpy.types.Collection.batoms = bpy.props.CollectionProperty(name = 'batoms', type = gui_blase.BlaseAtom)
    bpy.types.Collection.bond = bpy.props.CollectionProperty(name = 'bond', type = gui_blase.BlaseBond)
    bpy.types.Object.is_batom = bpy.props.BoolProperty(name = 'is_batom')
    bpy.types.Object.is_bcell = bpy.props.BoolProperty(name = 'is_bcell')
    bpy.types.Object.label = bpy.props.StringProperty(name = 'label')
    bpy.types.Object.species = bpy.props.StringProperty(name = 'species')
    bpy.types.Object.element = bpy.props.StringProperty(name = 'element')
    # bpy.types.Object.batom = bpy.props.PointerProperty(name = 'batom', type = gui_blase.BatomSettings)


    bpy.utils.register_class(gui_uilist.ListItem)
    bpy.utils.register_class(gui_uilist.MY_UL_List)
    bpy.utils.register_class(gui_uilist.LIST_OT_NewItem)
    bpy.utils.register_class(gui_uilist.LIST_OT_DeleteItem)
    bpy.utils.register_class(gui_uilist.LIST_OT_MoveItem)
    bpy.utils.register_class(gui_uilist.PT_ListExample)

    bpy.types.Scene.my_list = bpy.props.CollectionProperty(type = gui_uilist.ListItem)
    bpy.types.Scene.list_index = bpy.props.IntProperty(name = "Index for my_list",
                                             default = 0)



def unregister():
    bpy.types.TOPBAR_MT_file_import.remove(menu_func_import_blase)

    del bpy.types.Scene.my_list
    del bpy.types.Scene.list_index

    bpy.utils.unregister_class(gui_uilist.ListItem)
    bpy.utils.unregister_class(gui_uilist.MY_UL_List)
    bpy.utils.unregister_class(gui_uilist.LIST_OT_NewItem)
    bpy.utils.unregister_class(gui_uilist.LIST_OT_DeleteItem)
    bpy.utils.unregister_class(gui_uilist.LIST_OT_MoveItem)
    bpy.utils.unregister_class(gui_uilist.PT_ListExample)

    for cls in classes:
        bpy.utils.unregister_class(cls)


if __name__ == "__main__":

    register()
