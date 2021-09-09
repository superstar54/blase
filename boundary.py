# -*- coding: utf-8 -*-
"""
Author: 
    Xing Wang <xing.wang@psi.ch>
"""
import numpy as np
import math
import argparse
import time
import json
from ase.data import covalent_radii
from ase.visualize import view
from ase.neighborlist import NeighborList
from ase import Atom, Atoms
from blase.tools import euler_from_vector

class Boundary:
    """
    Boundary object.

    Example:
      cl = Boundary(atoms, d, index)
      cl.build()
    """

    def __init__(self, atoms, boundary_list, rotate_atoms = False):
        self.atoms = atoms
        self.natoms = len(atoms)
        self.cell = atoms.cell
        self.boundary_list = boundary_list
        # self.d = d
        # self.index = index
        self.rotate_atoms = rotate_atoms
        #
    def build(self, ):
        for boundary in self.boundary_list:
            print(boundary)
            self.cut(**boundary)
        return self.atoms
    def cut(self, atoms = None, d = None, index = None, direction = 1):
        """
        """
        if not atoms:
            atoms = self.atoms
        cell = atoms.cell
        normal = self.get_plane(d, index, cell)
        print(normal, d)
        # a*x + b*y + c*z - d = 0
        mask = []
        natoms = len(atoms)
        for i in range(natoms):
            v = atoms[i].position.dot(normal) - d
            # print(v)
            if v*direction > 0:
                mask.append(i)
        # view(atoms[mask])
        del atoms[mask]
        self.atoms = atoms
        if self.rotate_atoms:
            atoms = self.rotate()
    def rotate(self, atoms = None, index = None):
        """
        rotate normal of plane to z axis
        """
        import scipy
        if not atoms:
            atoms = self.atoms
        if not index:
            index = self.index
        cell = atoms.cell
        normal, d = self.get_plane(self.d, self.index, cell)
        vec = np.cross([0.0000014159, 0.000001951, 1], index)
        vec = vec/np.linalg.norm(vec)
        ang = np.arccos(normal[2])
        vec = ang*vec
        r = scipy.spatial.transform.Rotation.from_rotvec(vec)
        if scipy.version.version >= '1.4':
            mat = r.as_matrix()
        else: 
            mat = r.as_dcm()
        # print(mat)
        atoms.positions = atoms.positions.dot(mat)
        # atoms.cell = atoms.cell.dot(mat)
        self.atoms = atoms
        return atoms
    def get_plane(self, d, index = None, cell = None):
        '''
        plane equation: three point and distance from origin
        return normal and distance
        # a*x + b*y + c*z - d = 0
        '''
        index = [1.0/(index[i] + 0.000001) for i in range(3)]
        index = np.array(index)
        index = index/np.linalg.norm(index)
        points = cell*index
        # print(points)
        # x 
        v1 = points[1] - points[0]
        v2 = points[2] - points[0]
        normal = np.cross(v1, v2)
        print(v1, v2, normal)
        normal = normal/np.linalg.norm(normal)
        a, b, c = normal
        # d = np.dot(normal, points[2])
        return normal

if __name__ == "__main__":
    from ase.io import read, write
    from ase.build import molecule, bulk
    from ase.visualize import view
    from ase.data import covalent_radii
    # atoms = read('../examples/datas/tio2.cif')
    # atoms = read('../examples/datas/ceo2.cif')
    atoms = bulk('Pt', cubic = True)
    atoms = atoms*[6, 6, 6]
    # atoms.translate([0, 0, 0] - atoms.get_center_of_mass())
    # print(atoms)
    # cl = Boundary(atoms, d = 4.6,  index = [1, 1, 1], rotate_atoms = True)
    # cl = Boundary(atoms, boundary_list = [{'d': 10.0, 'index': [1, 1, 1], 'direction': -1},
                                          # {'d': 20.0, 'index': [1, 1, 1], 'direction': 1}]
        # )
    cl = Boundary(atoms, boundary_list = [{'d': 10.0, 'index': [2, 2, 1]}])

    # cl = Boundary(atoms, cutoffs = 1.2)
    atoms = cl.build()
    # atoms = cl.rotate()
    # print(atoms)
    view(atoms)

