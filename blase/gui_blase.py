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
                       BoolVectorProperty,
                       IntProperty,
                       IntVectorProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       EnumProperty,
                       PointerProperty,
                       )

from blase.gui_io import import_blase


from ase import Atom, Atoms
from ase.build import molecule, bulk
import json
from blase.bio import Blase
from blase.bdraw import draw_cell, draw_atoms, draw_bonds, draw_polyhedras, draw_isosurface, bond_source, cylinder_mesh_from_instance, clean_default
from blase.butils import read_blase_collection_list, read_blase_collection


class BlaseSettings(bpy.types.PropertyGroup):
    is_blase: BoolProperty(name="is_blase", default=False)
    model_type: StringProperty(name="model_type", default = '0')
    pbc: BoolVectorProperty(name="pbc", default = [False, False, False], size = 3)
    cell: FloatVectorProperty(name="cell", default = [0, 0, 0, 0, 0, 0, 0, 0, 0], size = 9)
    boundary: FloatVectorProperty(name="boundary", default = [0.0, 0.0, 0.0], size = 3)
class BlaseAtom(bpy.types.PropertyGroup):
    symbol: StringProperty(name="symbol")
    position: FloatVectorProperty(name="position", size = 3)
    tag: IntProperty(name="tag")
    
# The panel.
class Blase_PT_prepare(Panel):
    bl_label       = "Blase Tools"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    # bl_options     = {'DEFAULT_CLOSED'}
    bl_category = "Blase"
    bl_idname = "BLASE_PT_tools"

    filename: StringProperty(
        name = "Filename", default='blase-output.xyz',
        description = "Export atoms to file.")

    def draw(self, context):
        layout = self.layout
        blpanel = context.scene.blpanel

        box = layout.box()
        col = box.row()
        col.prop(blpanel, "collection_list")
        box = layout.box()
        col = box.column()
        col.label(text="Model")
        row = box.row()
        row.prop(blpanel, "model_type", expand  = True)
        box = layout.box()
        col = box.column()
        col.prop(blpanel, "boundary")

        box = layout.box()
        col = box.column(align=True)
        col.label(text="Atoms")
        row = box.row()
        row.prop(blpanel, "scale")

        

        box = layout.box()
        row = layout.row()
        row.prop(blpanel, "delete_index")

        box = layout.box()
        col = box.column(align=True)
        col.label(text="Split atoms")
        col.prop(blpanel, "single")
        row = box.row()
        row.operator("blase.split_atoms")

        box = layout.box()
        col = box.column(align=True)
        col.label(text="Copy atoms")
        col.prop(blpanel, "separate")
        row = box.row()
        row.operator("blase.copy_atoms")


        box = layout.box()
        col = box.column(align=True)
        col.label(text="Add structure")
        col.prop(blpanel, "atoms_str")
        col.prop(blpanel, "atoms_name")
        # col = box.column()
        col.prop(blpanel, "model_type_add")
        col.operator("blase.add_molecule")
        col.operator("blase.add_bulk")
        col.operator("blase.add_atoms")

        box = layout.box()
        col = box.column(align=True)
        col.label(text="Movie")
        col.prop(blpanel, "movie")
        
        box = layout.box()
        col = box.column(align=True)
        col.label(text="Render atoms")
        col.prop(blpanel, "output_image")



        box = layout.box()
        col = box.column(align=True)
        col.label(text="Export atoms")
        col.prop(blpanel, "output")
        col = box.column(align=True)
        col.operator("blase.export_atom")


class BlaseProperties(bpy.types.PropertyGroup):
    def Callback_model_type(self, context):
        blpanel = bpy.context.scene.blpanel
        print('Callback_model_type')
        model_type = list(blpanel.model_type)[0]
        modify_model_type(blpanel.collection_list, model_type)
    def Callback_delete_index(self, context):
        blpanel = bpy.context.scene.blpanel
        print('Callback_delete_index')
        modify_delete_index(blpanel.collection_list, blpanel.delete_index)
    def Callback_collection_list(self, context):
        print('Callback_collection_list')
        items = read_blase_collection_list()
        items = [(item, item, "") for item in items]
        items = tuple(items)
        return items
    def Callback_modify_scale(self, context):
        blpanel = bpy.context.scene.blpanel
        print('Callback_modify_scale')
        modify_scale(blpanel.collection_list, blpanel.scale)
    def Callback_modify_boundary(self, context):
        blpanel = bpy.context.scene.blpanel
        print('Callback_modify_boundary')
        modify_boundary(blpanel.collection_list, blpanel.boundary)
    
    def Callback_render_atoms(self, context):
        blpanel = bpy.context.scene.blpanel
        print('Callback_render_atoms')
        render_atoms(blpanel.collection_list, blpanel.output_image)
    def Callback_load_frames(self, context):
        blpanel = bpy.context.scene.blpanel
        print('Callback_load_frames')
        load_frames(blpanel.collection_list, blpanel.movie)


    
    collection_list: EnumProperty(
        name="Collection",
        description="Collection",
        items=Callback_collection_list)
    model_type: EnumProperty(
        name="Type",
        description="Structural models",
        items=(('0',"Space-filling", "Use ball and stick"),
               ('1',"Ball-and-stick", "Use ball"),
               ('2',"Polyhedral","Use polyhedral"),
               ('3',"Stick", "Use stick")),
        default={'0'}, 
        update=Callback_model_type,
        options={'ENUM_FLAG'},
        # options = {'ENUM_FLAG'}
        )
    model_type_add: EnumProperty(
        name="Type",
        description="Structural models",
        items=(('0',"Space-filling", "Use ball"),
               ('1',"Ball-and-stick", "Use ball and stick"),
               ('2',"Polyhedral","Use polyhedral"),
               ('3',"Stick", "Use stick")),
               default='0')
    delete_index: StringProperty(
        name="Delete",
        description="delete atoms",
        default='0', update=Callback_delete_index)
    action_type: EnumProperty(
        name="",
        description="Which objects shall be modified?",
        items=(('ALL_ACTIVE',"all active objects", "in the current layer"),
               ('ALL_IN_LAYER',"all in all selected layers",
                "in selected layer(s)")),
               default='ALL_ACTIVE',)
    scale: FloatProperty(
        name="scale", default=1.0,
        description = "scale", update = Callback_modify_scale)
    boundary: FloatVectorProperty(
        name="Boundary", default=(0.00, 0.0, 0.0),
        subtype = "XYZ",
        description = "boundary  in a, b, c axis", update = Callback_modify_boundary)
    single: BoolProperty(
        name="Single", default=False,
        description = "Do you split into single atoms?")
    separate: BoolProperty(
        name="Siparate", default=False,
        description = "Do you separate the copied atoms?")
    output: StringProperty(
        name = "Output", default='blase-output.xyz',
        description = "Output file")
    atoms_str: StringProperty(
        name = "Formula", default='H2O',
        description = "atoms_str")
    atoms_name: StringProperty(
        name = "Name", default='h2o',
        description = "name")
    output_image: StringProperty(
        name = "Output image", default='blase.png',
        description = "output render image", update = Callback_render_atoms)
    movie: IntProperty(
        name = "Load frames", default=1,
        description = "load frames", update = Callback_load_frames)


# Button for export atoms
class ExportAtom(Operator):
    bl_idname = "blase.export_atom"
    bl_label = "Save"
    bl_description = ("Save selected atoms to files")

    def execute(self, context):
        export_atom(context.scene.blpanel.output)
        return {'FINISHED'}

class AddMolecule(Operator):
    bl_idname = "blase.add_molecule"
    bl_label = "Add molecule"
    bl_description = ("Add molecule")
    def execute(self, context):
        #print(molecule_str)
        blpanel = context.scene.blpanel
        atoms = molecule(blpanel.atoms_str)
        import_blase(atoms, name = blpanel.atoms_name, model_type = blpanel.model_type_add)
        return {'FINISHED'}
class AddBulk(Operator):
    bl_idname = "blase.add_bulk"
    bl_label = "Add bulk"
    bl_description = ("Add bulk")
    def execute(self, context):
        #print(bulk_str)
        blpanel = context.scene.blpanel
        atoms = bulk(blpanel.atoms_str)
        import_blase(atoms, name = blpanel.atoms_name, model_type = blpanel.model_type_add, search_pbc_atoms={}, show_unit_cell=True)
        return {'FINISHED'}     
class AddAtoms(Operator):
    bl_idname = "blase.add_atoms"
    bl_label = "Add atoms"
    bl_description = ("Add atoms")
    def execute(self, context):
        #print(atoms_str)
        blpanel = context.scene.blpanel
        atoms = Atoms(blpanel.atoms_str)
        import_blase(atoms, name = blpanel.atoms_name, model_type = blpanel.model_type_add)
        return {'FINISHED'}     



class SplitAtom(Operator):
    bl_idname = "blase.split_atoms"
    bl_label = "Split"
    bl_description = ("Split atoms")

    def execute(self, context):
        split_atoms(context.scene.blpanel.single)
        return {'FINISHED'}
def split_atoms(single=False):
    '''
    '''
    objs = [obj for obj in bpy.context.selected_objects if obj.type == 'MESH']
    for obj in objs:
        mesh = obj.data
        if not obj.instance_type == "VERTS":
            return {'FINISHED'}
        coll_atom_kinds = obj.users_collection[0]
        obm = bmesh.new()
        if obj.mode == 'OBJECT':
            obm.from_mesh(mesh)
        elif obj.mode == 'EDIT':
            obm = bmesh.from_edit_mesh(mesh)
        positions = []
        for vert in obm.verts:
            if vert.select:
                positions.append(obj.matrix_world @ vert.co)
        obm.free()
        del(obm)
        sphere = obj.children[0]
        coll_instancers= sphere.users_collection[0]
        #
        print(coll_atom_kinds)
        if single:
            i = 1
            for position in positions:
                newsphere = sphere.copy()
                newsphere.data = sphere.data.copy()
                newsphere.name = '{0}_{1}'.format(sphere.name, i)
                # bpy.context.view_layer.objects.active = newsphere
                newmesh = bpy.data.meshes.new('{0}_{1}'.format(mesh.name, i))
                obj_atom = bpy.data.objects.new('{0}_{1}'.format(obj.name, i), newmesh)
                obj_atom.data.from_pydata([position], [], [])
                newsphere.parent = obj_atom
                obj_atom.instance_type = 'VERTS'
                coll_atom_kinds.objects.link(obj_atom)
                coll_instancers.objects.link(newsphere)
                newsphere.hide_set(True)
                i += 1
        else:
            newsphere = sphere.copy()
            newsphere.data = sphere.data.copy()
            newsphere.name = '{0}_1'.format(sphere.name)
            # bpy.context.view_layer.objects.active = newsphere
            newmesh = bpy.data.meshes.new('{0}_1'.format(mesh.name))
            obj_atom = bpy.data.objects.new('{0}_1'.format(obj.name), newmesh)
            obj_atom.data.from_pydata(positions, [], [])
            newsphere.parent = obj_atom
            obj_atom.instance_type = 'VERTS'
            coll_atom_kinds.objects.link(obj_atom)
            coll_instancers.objects.link(newsphere)
            newsphere.hide_set(True)
            # print(newsphere, obj_atom)
        if obj.mode == 'OBJECT':
            bpy.data.objects.remove(obj)
            bpy.data.objects.remove(sphere)
        elif obj.mode == 'EDIT':
            bpy.ops.mesh.delete(type='VERT')

class CopyAtoms(Operator):
    bl_idname = "blase.copy_atoms"
    bl_label = "Copy"
    bl_description = ("Copy atoms")

    def execute(self, context):
        copy_atoms(context.scene.blpanel.separate)
        return {'FINISHED'}

def copy_atoms(separate=False):
    '''
    '''
    #print(separate)
    # loop all selected object
    objs = [obj for obj in bpy.context.selected_objects if obj.type == 'MESH']
    for obj in objs:
        mesh = obj.data
        if not obj.instance_type == "VERTS":
            continue
        coll_atom_kinds = obj.users_collection[0]
        if obj.mode == 'OBJECT':
            obm = bmesh.new()
            obm.from_mesh(mesh)
        elif obj.mode == 'EDIT':
            obm = bmesh.from_edit_mesh(mesh)
        positions = []
        for vert in obm.verts:
            if vert.select:
                pos = obj.matrix_world @ vert.co + Vector([1, 1, 1])     
                positions.append(pos)
        
        sphere = obj.children[0]
        coll_instancers= sphere.users_collection[0]
        #
        print(coll_atom_kinds)
        print(positions)
        # merge the copied atoms to original atoms or not?
        if separate:
            # new atoms
            newsphere = sphere.copy()
            newsphere.data = sphere.data.copy()
            newsphere.name = '{0}_copy'.format(sphere.name)
            # bpy.context.view_layer.objects.active = newsphere
            newmesh = bpy.data.meshes.new('{0}_copy'.format(mesh.name))
            obj_atom = bpy.data.objects.new('{0}_copy'.format(obj.name), newmesh)
            obj_atom.data.from_pydata(positions, [], [])
            newsphere.parent = obj_atom
            obj_atom.instance_type = 'VERTS'
            coll_atom_kinds.objects.link(obj_atom)
            coll_instancers.objects.link(newsphere)
            obj_atom.select_set(True)
            obj.select_set(False)
            newsphere.hide_set(True)
        else:
            # merge
            # add to obm
            for pos in positions:
                obm.verts.new(pos)
            if obj.mode == 'EDIT':
                bmesh.update_edit_mesh(mesh)
                bpy.ops.object.mode_set(mode="OBJECT")
            else:
                # select the copied atoms, because we ar in "object" model, all vert are selected.
                for vert in obm.verts:
                    vert.select = True
                obm.to_mesh(mesh)
                mesh.update()        
            #print(len(mesh.vertices))
        #
        obm.free()
        del(obm)



def choose_objects(action_type,
                   model_type):

    # For selected objects of all selected layers
    change_objects = []
    # Note all selected objects first.
    for atom in bpy.context.selected_objects:
        change_objects.append(atom)
    print(change_objects)

    # This is very important now: If there are dupliverts structures, note
    # only the parents and NOT the children! Otherwise the double work is
    # done or the system can even crash if objects are deleted. - The
    # chidlren are accessed anyways (see below).
    change_objects = []
    for atom in change_objects_all:
        if atom.parent != None:
            FLAG = False
            for atom2 in change_objects:
                if atom2 == atom.parent:
                   FLAG = True
            if FLAG == False:
                change_objects.append(atom)
        else:
            change_objects.append(atom)

    # And now, consider all objects, which are in the list 'change_objects'.
    for atom in change_objects:
        if len(atom.children) != 0:
            for atom_child in atom.children:
                if atom_child.type in {'SURFACE', 'MESH', 'META'}:
                    modify_objects(action_type,
                                   atom_child,
                                   model_type)
        else:
            if atom.type in {'SURFACE', 'MESH', 'META'}:
                modify_objects(action_type,
                               atom,
                               model_type)



def read_atoms_select():
    '''   
    '''
    from ase import Atoms, Atom
    atoms = Atoms()
    cell_vertexs = []
    for obj in bpy.context.selected_objects:
        if "BOND" in obj.name.upper():
            continue
        if obj.type not in {'MESH', 'SURFACE', 'META'}:
            continue
        name = ""
        if 'atom_kind_' == obj.name[0:10]:
            print(obj.name)
            ind = obj.name.index('atom_kind_')
            ele = obj.name[ind + 10:].split('_')[0]
            if len(obj.children) != 0:
                for vertex in obj.data.vertices:
                    location = obj.matrix_world @ vertex.co
                    atoms.append(Atom(ele, location))
            else:
                if not obj.parent:
                    location = obj.location
                    atoms.append(Atom(ele, location))
        # cell
        if 'point_cell' == obj.name[0:10]:
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
        atoms.pbc = [True, True, True]
    # self.atoms = atoms
    return atoms
def export_atom(filename = 'test-blase.xyz'):
    atoms = read_atoms_select()
    atoms.write(filename)



# Modifying the scale of a selected atom or stick
def modify_model_type(collection_name, model_type, batoms = None):
    # Modify atom scale (all selected)
    coll = bpy.data.collections[collection_name]
    if not batoms:
        batoms = read_blase_collection(coll)
    for obj in coll.all_objects:
        bpy.data.objects.remove(obj)
    print('drawing atoms')
    batoms.coll.blase.model_type = model_type
    batoms.draw()
# Modifying the scale of a selected atom or stick
def modify_delete_index(collection_name, index, batoms = None):
    # Modify atom scale (all selected)
    coll = bpy.data.collections[collection_name]
    if not batoms:
        batoms = read_blase_collection(coll)
    print('Delete atoms')
    batoms.delete(index)
# Modifying the scale of a selected atom or stick
def modify_scale(collection_name, scale, batoms = None):
    # Modify atom scale (all selected)
    coll = bpy.data.collections[collection_name]
    for obj in coll.all_objects:
        if 'sphere' in obj.name:
            obj.scale = [scale, scale, scale]
# Modifying the scale of a selected atom or stick
def modify_boundary(collection_name, cutoff, batoms = None):
    # Modify atom cutoff (all selected)
    coll = bpy.data.collections[collection_name]
    if not batoms:
        batoms = read_blase_collection(coll)
    print('drawing atoms')
    batoms.coll.blase.boundary = cutoff
    # if batoms.atoms.pbc.any():
        # batoms.draw_atoms_boundary(cutoff=cutoff)
    batoms.draw()

def render_atoms(collection_name, output_image = 'bout.png', batoms = None):
    """
    """
    from blase.bio import Blase
    print('Rendering atoms')
    coll = bpy.data.collections[collection_name]
    if not batoms:
        batoms = read_blase_collection(coll)
    batoms.render(output_image = output_image)
def load_frames(collection_name, movie, batoms = None):
    """
    """
    from blase.bio import Blase
    print('Rendering atoms')
    coll = bpy.data.collections[collection_name]
    if not batoms:
        batoms = read_blase_collection(coll)
    batoms.load_frames()
