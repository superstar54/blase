import bpy
import numpy as np
from mathutils import Matrix
from scipy.spatial.transform import Rotation as R
import time
######################################################
#========================================================
def draw_cell(bobj = None, coll = None, cell_vertices = None, celllinewidth = None):
    """
    Draw unit cell
    """
    if not  coll:
        coll = bobj.coll
    coll_cell = [c for c in coll.children if c.name[0:4] == 'cell'][0]
    if not cell_vertices:
        cell_vertices = bobj.cell_vertices
    if not celllinewidth:
        celllinewidth = bobj.celllinewidth
    #------------------------------------------------------------
    if cell_vertices is not None:
        # build materials
        material = bpy.data.materials.new('cell')
        material.name = 'cell'
        material.diffuse_color = (0.8, 0.25, 0.25, 1.0)
        # draw points
        bpy.ops.mesh.primitive_uv_sphere_add(radius = celllinewidth) #, segments=32, ring_count=16)
        sphere = bpy.context.view_layer.objects.active
        sphere.name = 'sphere_cell'
        sphere.data.materials.append(material)
        bpy.ops.object.shade_smooth()
        sphere.hide_set(True)
        mesh = bpy.data.meshes.new('mesh_cell' )
        obj_cell = bpy.data.objects.new('point_cell', mesh )
        # Associate the vertices
        obj_cell.data.from_pydata(cell_vertices, [], [])
        # Make the object parent of the cube
        sphere.parent = obj_cell
        # Make the object dupliverts
        obj_cell.instance_type = 'VERTS'
        # STRUCTURE.append(obj_cell)
        bpy.data.collections['instancers'].objects.link(sphere)
        coll_cell.objects.link(obj_cell)
        #
        # edges
        edges = [[0, 1], [0, 2], [0, 4], 
                 [1, 3], [1, 5], [2, 3], 
                 [2, 6], [3, 7], [4, 5], 
                 [4, 6], [5, 7], [6, 7],
        ]
        cell_edges = {'lengths': [], 
                      'centers': [],
                      'normals': []}
        for e in edges:
            center = (cell_vertices[e[0]] + cell_vertices[e[1]])/2.0
            vec = cell_vertices[e[0]] - cell_vertices[e[1]]
            length = np.linalg.norm(vec)
            nvec = vec/length
            # print(center, nvec, length)
            cell_edges['lengths'].append(length/2.0)
            cell_edges['centers'].append(center)
            cell_edges['normals'].append(nvec)
        #
        source = bond_source(vertices=4)
        verts, faces = cylinder_mesh_from_instance(cell_edges['centers'], cell_edges['normals'], cell_edges['lengths'], celllinewidth, source)
        # print(verts)
        mesh = bpy.data.meshes.new("mesh_cell")
        mesh.from_pydata(verts, [], faces)  
        mesh.update()
        for f in mesh.polygons:
            f.use_smooth = True
        obj_edge = bpy.data.objects.new("edge_cell", mesh)
        obj_edge.data = mesh
        obj_edge.data.materials.append(material)
        bpy.ops.object.shade_smooth()
        coll_cell.objects.link(obj_edge)

def draw_atoms(bobj = None, coll = None, atom_kinds = None, bsdf_inputs = None, material_style = 'blase', make_real = None):
    '''
    Draw atoms
    bsdf_inputs: dict
        The key and value for principled_bsdf node
    material_style: string
        Select materials type from ['blase', 'glass', 'ceramic', 'plastic'].
    '''
    # build materials
    if not  coll:
        coll = bobj.coll
    coll_atom_kinds = [c for c in coll.children if 'atoms' in c.name][0]
    if not atom_kinds:
        atom_kinds = bobj.atom_kinds
    if not bsdf_inputs:
        bsdf_inputs = bobj.material_styles_dict[material_style]
    if make_real is None:
        make_real = bobj.make_real
    for kind, datas in atom_kinds.items():
        tstart = time.time()
        material = bpy.data.materials.new('atom_kind_{0}'.format(kind))
        material.diffuse_color = np.append(datas['color'], datas['transmit'])
        # material.blend_method = 'BLEND'
        material.use_nodes = True
        principled_node = material.node_tree.nodes['Principled BSDF']
        principled_node.inputs['Base Color'].default_value = np.append(datas['color'], datas['transmit'])
        principled_node.inputs['Alpha'].default_value = datas['transmit']
        for key, value in bsdf_inputs.items():
            principled_node.inputs[key].default_value = value
        datas['materials'] = material
        #
        # print(datas['positions'][0])
        if datas['balltype'] == 'meta':
            print('metaball', kind)
            mbdata = bpy.data.metaballs.new('atom_kind'.format(kind))
            mbdata.render_resolution = 0.1
            mbdata.resolution = 0.2
            obj_atom = bpy.data.objects.new('meta_atom_kind'.format(kind), mbdata)
            obj_atom.data.materials.append(material)
            for co in datas['positions']:
                mbele = mbdata.elements.new(type = 'BALL')
                mbele.co = co
                mbele.radius = datas['radius']*2.0
                mbele.stiffness = 1.0
            # STRUCTURE.append(obj_atom)
            coll_atom_kinds.objects.link(obj_atom)
        else:
            bpy.ops.mesh.primitive_uv_sphere_add(radius = datas['radius']) #, segments=32, ring_count=16)
            sphere = bpy.context.view_layer.objects.active
            sphere.name = 'sphere_atom_kind_{0}'.format(kind)
            sphere.data.materials.append(material)
            bpy.ops.object.shade_smooth()
            sphere.hide_set(True)
            mesh = bpy.data.meshes.new('mesh_kind_{0}'.format(kind) )
            obj_atom = bpy.data.objects.new('atom_kind_{0}'.format(kind), mesh )
            # Associate the vertices
            obj_atom.data.from_pydata(datas['positions'], [], [])
            # Make the object parent of the cube
            sphere.parent = obj_atom
            # Make the object dupliverts
            obj_atom.instance_type = 'VERTS'
            # bpy.context.view_layer.objects.active = obj_atom
            # obj_atom.select_set(True)
            # bpy.data.objects['atom_kind_{0}'.format(kind)].select_set(True)
            # STRUCTURE.append(obj_atom)
            bpy.data.collections['instancers'].objects.link(sphere)
            coll_atom_kinds.objects.link(obj_atom)
        print('atoms: {0}   {1:10.2f} s'.format(kind, time.time() - tstart))
    if make_real:
        tstart = time.time()
        bpy.ops.object.select_by_type(extend=False, type='MESH')
        bpy.ops.object.duplicates_make_real()
        for kind in atom_kinds:
            bpy.data.objects.remove(bpy.data.objects['atom_kind_{0}'.format(kind)])
        print('make_real: {0:10.2f} s'.format(time.time() - tstart))

def draw_bonds(bobj = None, coll = None, bond_kinds = None, bond_list= None, bondlinewidth = None, vertices = None, bsdf_inputs = None, material_style = 'blase'):
    '''
    Draw atom bonds
    '''
    if not  coll:
        coll = bobj.coll
    coll_bond_kinds = [c for c in coll.children if 'bonds' in c.name][0]
    if not bond_kinds:
        bond_kinds = bobj.bond_kinds
    if not bondlinewidth:
        bondlinewidth = bobj.bondlinewidth
    if not bsdf_inputs:
        bsdf_inputs = bobj.material_styles_dict[material_style]
    # import pprint
    # pprint.pprint(bond_kinds)
    if not vertices:
    	if len(bobj.atoms) > 200:
	        vertices=8
    	else:
    		vertices = 16
    source = bond_source(vertices = vertices)
    #
    for kind, datas in bond_kinds.items():
        tstart = time.time()
        material = bpy.data.materials.new('bond_kind_{0}'.format(kind))
        material.diffuse_color = np.append(bond_kinds[kind]['color'], bond_kinds[kind]['transmit'])
        # material.blend_method = 'BLEND'
        material.use_nodes = True
        principled_node = material.node_tree.nodes['Principled BSDF']
        principled_node.inputs['Base Color'].default_value = np.append(datas['color'], datas['transmit'])
        principled_node.inputs['Alpha'].default_value = bond_kinds[kind]['transmit']
        for key, value in bsdf_inputs.items():
            principled_node.inputs[key].default_value = value
        datas['materials'] = material
        #
        # print(datas['normals'])
        verts, faces = cylinder_mesh_from_instance(datas['centers'], datas['normals'], datas['lengths'], bondlinewidth, source)
        mesh = bpy.data.meshes.new("mesh_kind_{0}".format(kind))
        mesh.from_pydata(verts, [], faces)  
        mesh.update()
        for f in mesh.polygons:
            f.use_smooth = True
        obj_bond = bpy.data.objects.new("bond_kind_{0}".format(kind), mesh)
        obj_bond.data = mesh
        obj_bond.data.materials.append(material)
        bpy.ops.object.shade_smooth()
        coll_bond_kinds.objects.link(obj_bond)
        print('bonds: {0}   {1:10.2f} s'.format(kind, time.time() - tstart))

def draw_bonds_2(bobj = None, coll = None, bond_kinds = None, bond_list= None, bondlinewidth = None, vertices = None, bsdf_inputs = None, material_style = 'blase'):
    '''
    Draw atom bonds
    '''
    if not  coll:
        coll = bobj.coll
    coll_bond_kinds = [c for c in coll.children if 'bonds' in c.name][0]
    if not bond_kinds:
        bond_kinds = bobj.bond_kinds
    if not bondlinewidth:
        bondlinewidth = bobj.bondlinewidth
    if not bsdf_inputs:
        bsdf_inputs = bobj.material_styles_dict[material_style]
    # import pprint
    # pprint.pprint(bond_kinds)
    if not vertices:
        if len(bobj.atoms) > 200:
            vertices=8
        else:
            vertices = 16
    #
    # print(bond_kinds)
    for kind, datas in bond_kinds.items():
        tstart = time.time()
        material = bpy.data.materials.new('bond_kind_{0}'.format(kind))
        material.diffuse_color = np.append(bond_kinds[kind]['color'], bond_kinds[kind]['transmit'])
        # material.blend_method = 'BLEND'
        material.use_nodes = True
        principled_node = material.node_tree.nodes['Principled BSDF']
        principled_node.inputs['Base Color'].default_value = np.append(datas['color'], datas['transmit'])
        principled_node.inputs['Alpha'].default_value = bond_kinds[kind]['transmit']
        for key, value in bsdf_inputs.items():
            principled_node.inputs[key].default_value = value
        datas['materials'] = material
        #
        bpy.ops.mesh.primitive_cylinder_add(vertices = vertices, radius=0.1, depth = 1.0)
        cylinder = bpy.context.view_layer.objects.active
        cylinder.name = 'cylinder_atom_kind_{0}'.format(kind)
        cylinder.data.materials.append(material)
        bpy.ops.object.shade_smooth()
        cylinder.hide_set(True)
        mesh = bpy.data.meshes.new('mesh_kind_{0}'.format(kind) )
        obj_bond = bpy.data.objects.new('bond_kind_{0}'.format(kind), mesh )
        # Associate the vertices
        obj_bond.data.from_pydata(datas['verts'], [], datas['faces'])
        # Make the object parent of the cube
        cylinder.parent = obj_bond
        # Make the object dupliverts
        obj_bond.instance_type = 'FACES'
        obj_bond.use_instance_faces_scale = True
        obj_bond.show_instancer_for_render = False
        obj_bond.show_instancer_for_viewport = False
        # bpy.context.view_layer.objects.active = obj_bond
        # obj_bond.select_set(True)
        # bpy.data.objects['bond_kind_{0}'.format(kind)].select_set(True)
        # STRUCTURE.append(obj_bond)
        bpy.data.collections['instancers'].objects.link(cylinder)
        coll_bond_kinds.objects.link(obj_bond)        
        print('bonds: {0}   {1:10.2f} s'.format(kind, time.time() - tstart))
def draw_polyhedras(bobj, coll = None, polyhedra_kinds = None, polyhedra_dict= None, bsdf_inputs = None, material_style = 'plastic'):
    '''
    Draw polyhedras
    '''
    if not  coll:
        coll = bobj.coll
    coll_polyhedra_kinds = [c for c in coll.children if 'polyhedras' in c.name][0]
    if not polyhedra_kinds:
        polyhedra_kinds = bobj.polyhedra_kinds
    if not bsdf_inputs:
        bsdf_inputs = bobj.material_styles_dict[material_style]
    #
    
    # import pprint
    # pprint.pprint(polyhedra_kinds)
    source = bond_source(vertices=4)
    for kind, datas in polyhedra_kinds.items():
        tstart = time.time()
        material = bpy.data.materials.new('polyhedra_kind_{0}'.format(kind))
        material.diffuse_color = np.append(datas['color'], datas['transmit'])
        # material.blend_method = 'BLEND'
        material.use_nodes = True
        principled_node = material.node_tree.nodes['Principled BSDF']
        principled_node.inputs['Base Color'].default_value = np.append(datas['color'], datas['transmit'])
        principled_node.inputs['Alpha'].default_value = polyhedra_kinds[kind]['transmit']
        for key, value in bsdf_inputs.items():
            principled_node.inputs[key].default_value = value
        datas['materials'] = material
        #
        # create new mesh structure
        mesh = bpy.data.meshes.new("mesh_kind_{0}".format(kind))
        # mesh.from_pydata(polyhedra_kinds[kind]['vertices'], polyhedra_kinds[kind]['edges'], polyhedra_kinds[kind]['faces'])  
        mesh.from_pydata(datas['vertices'], [], datas['faces'])  
        mesh.update()
        for f in mesh.polygons:
            f.use_smooth = True
        obj_polyhedra = bpy.data.objects.new("polyhedra_kind_{0}".format(kind), mesh)
        obj_polyhedra.data = mesh
        obj_polyhedra.data.materials.append(material)
        bpy.ops.object.shade_smooth()
        #---------------------------------------------------
        material = bpy.data.materials.new('polyhedra_edge_kind_{0}'.format(kind))
        material.diffuse_color = np.append(datas['edge_cylinder']['color'], datas['edge_cylinder']['transmit'])
        # material.blend_method = 'BLEND'
        material.use_nodes = True
        principled_node = material.node_tree.nodes['Principled BSDF']
        principled_node.inputs['Base Color'].default_value = np.append(datas['edge_cylinder']['color'], datas['edge_cylinder']['transmit'])
        principled_node.inputs['Alpha'].default_value = datas['transmit']
        for key, value in bsdf_inputs.items():
            principled_node.inputs[key].default_value = value
        datas['edge_cylinder']['materials'] = material
        verts, faces = cylinder_mesh_from_instance(datas['edge_cylinder']['centers'], datas['edge_cylinder']['normals'], datas['edge_cylinder']['lengths'], 0.01, source)
        # print(verts)
        mesh = bpy.data.meshes.new("mesh_kind_{0}".format(kind))
        mesh.from_pydata(verts, [], faces)  
        mesh.update()
        for f in mesh.polygons:
            f.use_smooth = True
        obj_edge = bpy.data.objects.new("edge_kind_{0}".format(kind), mesh)
        obj_edge.data = mesh
        obj_edge.data.materials.append(material)
        bpy.ops.object.shade_smooth()
        # STRUCTURE.append(obj_polyhedra)
        coll_polyhedra_kinds.objects.link(obj_polyhedra)
        coll_polyhedra_kinds.objects.link(obj_edge)
        print('polyhedras: {0}   {1:10.2f} s'.format(kind, time.time() - tstart))

def draw_isosurface(bobj = None, coll = None, volume = None, level = 0.02,
                    closed_edges = False, gradient_direction = 'descent',
                    color=(0.85, 0.80, 0.25) , icolor = None, transmit=0.5,
                    verbose = False, step_size = 1, 
                    bsdf_inputs = None, material_style = 'blase'):
    """Computes an isosurface from a volume grid.
    
    Parameters:     
    """
    from skimage import measure
    colors = [(0.85, 0.80, 0.25), (0.0, 0.0, 1.0)]
    if icolor:
        color = colors[icolor]
    if not  coll:
        coll = bobj.coll
    coll_isosurface = [c for c in coll.children if 'isosurfaces' in c.name][0]
    
    cell = bobj.cell
    bobj.cell_vertices.shape = (2, 2, 2, 3)
    cell_origin = bobj.cell_vertices[0,0,0]
    #
    spacing = tuple(1.0/np.array(volume.shape))
    scaled_verts, faces, normals, values = measure.marching_cubes_lewiner(volume, level = level,
                    spacing=spacing,gradient_direction=gradient_direction , 
                    allow_degenerate = False, step_size=step_size)
    #
    scaled_verts = list(scaled_verts)
    nverts = len(scaled_verts)
    # transform
    for i in range(nverts):
        scaled_verts[i] = scaled_verts[i].dot(cell)
        scaled_verts[i] -= cell_origin
    faces = list(faces)
    print('Draw isosurface...')
    # print('verts: ', scaled_verts[0:5])
    # print('faces: ', faces[0:5])
    #material
    if not bsdf_inputs:
        bsdf_inputs = bobj.material_styles_dict[material_style]
    material = bpy.data.materials.new('isosurface')
    material.name = 'isosurface'
    material.diffuse_color = color + (transmit,)
    # material.alpha_threshold = 0.2
    # material.blend_method = 'BLEND'
    material.use_nodes = True
    principled_node = material.node_tree.nodes['Principled BSDF']
    principled_node.inputs['Base Color'].default_value = color + (transmit,)
    principled_node.inputs['Alpha'].default_value = transmit
    for key, value in bsdf_inputs.items():
            principled_node.inputs[key].default_value = value
    #
    # create new mesh structure
    isosurface = bpy.data.meshes.new("isosurface")
    isosurface.from_pydata(scaled_verts, [], faces)  
    isosurface.update()
    for f in isosurface.polygons:
        f.use_smooth = True
    iso_object = bpy.data.objects.new("isosurface", isosurface)
    iso_object.data = isosurface
    iso_object.data.materials.append(material)
    bpy.ops.object.shade_smooth()
    coll_isosurface.objects.link(iso_object)

def clean_default():
    if 'Camera' in bpy.data.cameras:
        bpy.data.cameras.remove(bpy.data.cameras['Camera'])
    if 'Light' in bpy.data.lights:
        bpy.data.lights.remove(bpy.data.lights['Light'])
    if 'Cube' in bpy.data.objects:
        bpy.data.objects.remove(bpy.data.objects['Cube'])



# draw bonds
def bond_source(vertices = 16):
    bpy.ops.mesh.primitive_cylinder_add(vertices = vertices)
    cyli = bpy.context.view_layer.objects.active
    me = cyli.data
    verts = []
    faces = []
    for vertices in me.vertices:
        verts.append(np.array(vertices.co))
    for poly in me.polygons:
        face = []
        for loop_index in range(poly.loop_start, poly.loop_start + poly.loop_total):
            # print("    Vertex: %d" % me.loops[loop_index].vertex_index)
            face.append(me.loops[loop_index].vertex_index)
        faces.append(face)
    cyli.select_set(True)
    bpy.ops.object.delete()
    return [verts, faces]
# draw atoms
def atom_source():
    bpy.ops.mesh.primitive_uv_sphere_add() #, segments=32, ring_count=16)
    # bpy.ops.mesh.primitive_cylinder_add()
    sphe = bpy.context.view_layer.objects.active
    me = sphe.data
    verts = []
    faces = []
    for vertices in me.vertices:
        verts.append(np.array(vertices.co))
    for poly in me.polygons:
        face = []
        for loop_index in range(poly.loop_start, poly.loop_start + poly.loop_total):
            # print("    Vertex: %d" % me.loops[loop_index].vertex_index)
            face.append(me.loops[loop_index].vertex_index)
        faces.append(face)
    sphe.select_set(True)
    bpy.ops.object.delete()
    return [verts, faces]



def sphere_mesh_from_instance(centers, radius, source):
    verts = []
    faces = []
    vert0, face0 = source
    nvert = len(vert0)
    nb = len(centers)
    for i in range(nb):
        center = centers[i]
        # r = radius[i]
        # normal = normal/np.linalg.norm(normal)
        for vert in vert0:
            vert = vert*[radius, radius, radius]
            vert += center
            verts.append(vert)
        for face in face0:
            face = [x+ i*nvert for x in face]
            faces.append(face)
    return verts, faces
def cylinder_mesh_from_instance(centers, normals, lengths, scale, source):
    # verts = np.empty((0, 3), float)
    verts = []
    faces = []
    vert0, face0 = source
    nvert = len(vert0)
    nb = len(centers)
    for i in range(nb):
        center = centers[i]
        normal = normals[i]
        length = lengths[i]
        # normal = normal/np.linalg.norm(normal)
        vec = np.cross([0.0000014159, 0.000001951, 1], normal)
        # print(center, normal, vec)
        vec = vec/np.linalg.norm(vec)
        # print(vec)
        # ang = np.arcsin(np.linalg.norm(vec))
        ang = np.arccos(normal[2])
        vec = -1*ang*vec
        r = R.from_rotvec(vec)
        matrix = r.as_dcm()
        # print(vec, ang)
        vert1 = vert0.copy()
        vert1 = vert1*np.array([scale, scale, length])
        vert1 = vert1.dot(matrix)
        vert1 += center
        for vert in vert1:
            verts.append(vert)
        for face in face0:
            face = [x+ i*nvert for x in face]
            faces.append(face)
    return verts, faces

