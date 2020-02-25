import bpy
import numpy as np
from mathutils import Matrix
from scipy.spatial.transform import Rotation as R
######################################################

def clean_default():
    bpy.data.cameras.remove(bpy.data.cameras['Camera'])
    bpy.data.lights.remove(bpy.data.lights['Light'])
    bpy.data.objects.remove(bpy.data.objects['Cube'])


def clean_objects():
    for item in bpy.data.objects:
        bpy.data.objects.remove(item)

def removeAll():
    #types =  ['MESH', 'CURVE', 'SURFACE', 'META', 'FONT', 'ARMATURE', 'LATTICE', 'EMPTY', 'GPENCIL', 'CAMERA', 'LIGHT', 'SPEAKER', 'LIGHT_PROBE']
    for mesh in bpy.data.meshes:
        bpy.data.meshes.remove(mesh)
    for obj in bpy.data.objects:
        bpy.data.objects.remove(obj)
    for cam in bpy.data.cameras:
        bpy.data.cameras.remove(cam)
    for light in bpy.data.lights:
        bpy.data.lights.remove(light)
    for coll in bpy.data.collections:
        bpy.data.collections.remove(coll)

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
        for vert in vert0:
            vert = vert*[scale, scale, length]
            vert = vert.dot(matrix)
            vert += center
            verts.append(vert)
        for face in face0:
            face = [x+ i*nvert for x in face]
            faces.append(face)
    return verts, faces

