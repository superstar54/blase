import bpy
from bpy.props import (StringProperty,
                       BoolProperty,
                       BoolVectorProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       EnumProperty,
                       PointerProperty,
                       )

class BlaseBatoms(bpy.types.PropertyGroup):
    is_batoms: BoolProperty(name="is_batoms", default=False)
    model_type: StringProperty(name="model_type", default = '0')
    pbc: BoolVectorProperty(name="pbc", default = [False, False, False], size = 3)
    cell: FloatVectorProperty(name="cell", default = [0, 0, 0, 0, 0, 0, 0, 0, 0], size = 9)
    show_unit_cell: BoolProperty(name="show_unit_cell", default = True)
    boundary: FloatVectorProperty(name="boundary", default = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0], size = 6)

class BlaseBatom(bpy.types.PropertyGroup):
    is_batom: BoolProperty(name="is_batom", default=False)
    label: StringProperty(name="label", default = '')
    species: StringProperty(name="species", default = 'X')
    element: StringProperty(name="element", default = '')
class BlaseBcell(bpy.types.PropertyGroup):
    is_bcell: BoolProperty(name="is_bcell", default=False)
    label: StringProperty(name="label", default = '')
class BlaseAtom(bpy.types.PropertyGroup):
    symbol: StringProperty(name="symbol")
    position: FloatVectorProperty(name="position", size = 3)
    tag: IntProperty(name="tag")
class BlaseBond(bpy.types.PropertyGroup):
    symbol1: StringProperty(name="symbol1")
    symbol2: StringProperty(name="symbol2")
    name:StringProperty(name = "name")
    min: FloatProperty(name="min", description = "min", default = 0.0)
    max: FloatProperty(name="max", description = "max", default = 2.0)
    polyhedra: BoolProperty(name="polyhedra", default=False)
    search: IntProperty(name="search", default=0)
    bondcolor1: FloatVectorProperty(name="bondcolor1", size = 4)
    bondcolor2: FloatVectorProperty(name="bondcolor1", size = 4)
    color2: FloatVectorProperty(name="bondcolor1", size = 4)
    bondlinewidth: FloatProperty(name="bondlinewidth", default = 0.10)
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Bondpair      min     max   Search_bond    Polyhedra \n'
        s += '{0:10s} {1:4.3f}   {2:4.3f}      {3:10s}   {4:10s} \n'.format(\
                self.name, self.min, self.max, str(self.search), str(self.polyhedra))
        s += '-'*60 + '\n'
        return s
