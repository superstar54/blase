# Blase TOOLBAR  - Addon - Blender 2.9x



###
bl_info = {
    "name": "Blase toolbar",
    "author": "Xing Wang",
    "version": (0, 5),
    "blender": (2, 93, 0),
    "location": "File -> Import -> Blase (xyz, cif, pdb, ...)",
    "description": "Python module for drawing and rendering ASE (Atomic Simulation Environment) atoms and molecules objects using blender.",
    "warning": "",
    "category": "Import-Export",
}


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
#
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
    

def unregister():
    bpy.types.TOPBAR_MT_file_import.remove(menu_func_import_blase)


    for cls in classes:
        bpy.utils.unregister_class(cls)


if __name__ == "__main__":

    register()
