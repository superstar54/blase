import bpy
from blase.batoms import Batoms
from ase import Atom, Atoms

def read_blase_collection_list():
    items = [col.name for col in bpy.data.collections if col.blase.is_blase]
    return items
def load_atoms(names = None):
    if names is None:
        names = read_blase_collection_list()
    objs = {}
    for name in names:
        coll = bpy.data.collections[name]
        atoms = read_blase_collection(coll)
        objs[name] = atoms
    return objs
def read_atoms_list(coll):
    '''   
    '''
    coll_atom_kinds = [coll for coll in coll.children if 'atoms' in coll.name][0]
    elements = []
    for obj in coll_atom_kinds.all_objects:
        ele = obj.name.split('_')[2]
        elements.append(ele)
    return elements
        
def read_blase_collection(coll):
    '''   
    '''
    atoms = Atoms()
    name = coll.name

    # atoms
    for obj in coll.children['%s_atoms'%name].all_objects:
        ele = obj.name.split('_')[2]
        if len(obj.children) != 0:
            for vertex in obj.data.vertices:
                location = obj.matrix_world @ vertex.co
                atoms.append(Atom(symbol = ele, position = location))
        else:
            if not obj.parent:
                location = obj.location
                atoms.append(Atom(ele, location))
    # atoms properties
    scale = {}
    for obj in coll.children['%s_instancers'%name].all_objects:
        ele = obj.name.split('_')[3]
        scale[ele] = obj.scale
    # cell
    coll_cell = coll.children['%s_cell'%name]
    cell_vertexs = []
    if 'point_cell' in coll_cell.all_objects.keys():
        obj = coll_cell.all_objects['point_cell']
        for vertex in obj.data.vertices:
            location = obj.matrix_world @ vertex.co
            cell_vertexs.append(location)
    if cell_vertexs:
        cell = [cell_vertexs[4], cell_vertexs[2], cell_vertexs[1]]
        atoms.cell = cell
        atoms.pbc = [True, True, True]
    # self.atoms = atoms
    batoms = Batoms(atoms, name = name, coll = coll, scale = scale)
    return batoms
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