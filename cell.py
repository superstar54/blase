"""
"""
from numpy.core.records import array
import bpy
import numpy as np
from ase.cell import Cell
from blase.btools import object_mode

class Bcell():
    """
    
    """
    def __init__(self, label, array = np.zeros([3, 3]), location = np.array([0, 0, 0])) -> None:
        """
        ver: 3x3 verlike object
          The three cell vectors: cell[0], cell[1], and cell[2].
        """
        self.label = label
        self.name = 'cell_%s_edge'%(self.label)
        self.draw_cell_edge(array, location)
    def draw_cell_edge(self, array, location):
        """
        Draw unit cell by edge, however, can not be rendered.
        """
        if array is None:
            array = np.zeros([3, 3])
        if not len(array) == 8:
            cell = Cell.new(array)
            verts = self.array2verts(cell.array)
        else:
            verts = array - array[3]
            location = array[3]
        if self.name not in bpy.data.objects:
            edges = [[3, 0], [3, 1], [4, 0], [4, 1],
                    [2, 5], [2, 6], [7, 5], [7, 6], 
                    [3, 2], [0, 6], [1, 5], [4, 7]
            ]
            mesh = bpy.data.meshes.new("cell_%s_edge"%self.label)
            mesh.from_pydata(verts, edges, [])  
            mesh.update()
            for f in mesh.polygons:
                f.use_smooth = True
            obj_edge = bpy.data.objects.new("cell_%s_edge"%self.label, mesh)
            obj_edge.data = mesh
            obj_edge.location = location
            bpy.data.collections['Collection'].objects.link(obj_edge)
        elif hasattr(bpy.data.objects[self.name], 'bcell'):
            obj_atom = bpy.data.objects[self.name]
        else:
            raise Exception("Failed, the name %s already in use and is not Bcell object!"%self.name)
        bpy.context.view_layer.update()
    def __repr__(self) -> str:
        numbers = self.array.tolist()
        s = 'Cell({})'.format(numbers)
        return s
    def __getitem__(self, index):
        return self.array[index]
    def __setitem__(self, index, value):
        """
        Add bondpair one by one
        """
        """Set unit cell vectors.

        Parameters:

        Examples:

        """
        from ase.cell import Cell
        from ase.geometry.cell import complete_cell
        bcell = self.bcell
        array = self.array
        array[index] = value
        verts = self.array2verts(array)
        for i in range(8):
            bcell.data.vertices[i].co = np.array(verts[i])
    def __array__(self, dtype=float):
        if dtype != float:
            raise ValueError('Cannot convert cell to array of type {}'
                             .format(dtype))
        return self.array
    @property
    def bcell(self):
        return self.get_bcell()
    def get_bcell(self):
        return bpy.data.objects['cell_%s_edge'%(self.label)]
    @property
    def array(self):
        return self.get_array()
    def get_array(self):
        cell = np.array([self.local_verts[0] - self.local_verts[3],
                         self.local_verts[1] - self.local_verts[3],
                         self.local_verts[2] - self.local_verts[3]])
        return cell
    @property
    def local_verts(self):
        return self.get_local_verts()
    def get_local_verts(self):
        bcell = self.bcell
        return np.array([bcell.data.vertices[i].co for i in range(8)])
    @property
    def verts(self):
        return self.get_verts()
    def get_verts(self):
        return np.array([self.bcell.matrix_world @ self.bcell.data.vertices[i].co for i in range(8)])
    def array2verts(self, array):
        """
        """
        verts = np.array([[1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [0, 0, 0],
            [1, 1, 0],
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 1],
            ])
        verts = np.dot(verts, array)
        return verts
    @property
    def location(self):
        return self.get_location()
    def get_location(self):
        return np.array(self.bcell.location)
    def copy(self, label):
        object_mode()
        cell = Bcell(label, array = self.array, location = self.bcell.location)
        return cell
    def repeat(self, m):
        self[:] = np.array([m[c] * self.array[c] for c in range(3)])
    def draw_cell_cylinder(self, celllinewidth = 0.01):
        """
        Draw unit cell using cylinder.
        """
        if self.local_verts is not None:
            # build materials
            material = bpy.data.materials.new('cell_{0}'.format(label))
            # material.label = 'cell'
            material.diffuse_color = (0.8, 0.25, 0.25, 1.0)
            # draw points
            bpy.ops.mesh.primitive_uv_sphere_add(radius = celllinewidth) #, segments=32, ring_count=16)
            sphere = bpy.context.view_layer.objects.active
            sphere.name = 'instancer_cell_%s_sphere'%label
            sphere.data.materials.append(material)
            bpy.ops.object.shade_smooth()
            sphere.hide_set(True)
            mesh = bpy.data.meshes.new('point_cell' )
            obj_cell = bpy.data.objects.new('cell_%s_point'%label, mesh )
            # Associate the vertices
            obj_cell.data.from_pydata(self.local_verts, [], [])
            sphere.parent = obj_cell
            obj_cell.instance_type = 'VERTS'
            coll_cell.objects.link(sphere)
            coll_cell.objects.link(obj_cell)
            #
            # edges
            edges = [[3, 0], [3, 1], [4, 0], [4, 1],
                    [2, 5], [2, 6], [7, 5], [7, 6], 
                    [3, 2], [0, 6], [1, 5], [4, 7]
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
            mesh = bpy.data.meshes.new("edge_cell")
            mesh.from_pydata(verts, [], faces)  
            mesh.update()
            for f in mesh.polygons:
                f.use_smooth = True
            obj_edge = bpy.data.objects.new("cell_%s_edge"%label, mesh)
            obj_edge.data = mesh
            obj_edge.data.materials.append(material)
            bpy.ops.object.shade_smooth()
            coll_cell.objects.link(obj_edge)
    