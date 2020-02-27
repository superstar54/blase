# -*- coding: utf-8 -*-
"""
This module solves the folllowing probllem:
Given the coordinates of a set of points on the plane, consider a link between 
two points if their distance is below a given threshold (obtaining in this way
an undirected graph). Then, given a list of pairs of points in input, return 
whether each pair of points is connected or not through a path in the graph. 

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

class CutSurface:
    """
    Connectivity list object.

    threshold: float
        If the distance of two points smaller than the threshold, they will be counted as
        connected. The (Cartesian) distance is defined as: d(p, q) = sqrt((px − qx)**2 + (py − qy)**2).
    

    Example:
      cl = CutSurface(threshold)
      cl.build(points)
      output = cl.get_connectivity(pairs)
    """

    def __init__(self, atoms, d, index = None):
        self.atoms = atoms
        self.natoms = len(atoms)
        self.cell = atoms.cell
        self.d = d
        self.index = index
        #
    def build(self, atoms = None, index = None):
        """
        """
        if not atoms:
            atoms = self.atoms
        natoms = len(atoms)
        d=10
        cell = atoms.cell
        index = np.array(self.index)
        index = index/np.linalg.norm(index)
        self.index = index
        points = cell*index*self.d
        print(points)
        # x 
        v1 = points[1] - points[0]
        v2 = points[2] - points[0]
        cp = np.cross(v1, v2)
        a, b, c = cp
        d = np.dot(cp, points[2])
        # a*x + b*y + c*z - d = 0
        mask = []
        for i in range(natoms):
            v = atoms[i].position.dot(cp) - d
            if v > 0:
                mask.append(i)
        # view(atoms[mask])
        del atoms[mask]
        self.atoms = atoms
        return atoms
    def rotate(self, atoms = None, index = None):
        """
        """
        if not atoms:
            atoms = self.atoms
        if not index:
            index = self.index
        index = np.array(self.index)
        ang = euler_from_vector(index, 'zyz')
        ang = ang*180/np.pi
        print(ang)
        atoms.euler_rotate(ang[0], ang[1], ang[2],)
        self.atoms = atoms
        return atoms

if __name__ == "__main__":
    from ase.io import read, write
    from ase.build import molecule, bulk
    from ase.visualize import view
    from ase.data import covalent_radii
    # atoms = read('../examples/datas/tio2.cif')
    # atoms = read('../examples/datas/ceo2.cif')
    atoms = bulk('Pt', cubic = True)
    atoms = atoms*[4, 4, 4]
    # print(atoms)
    cl = CutSurface(atoms, d = 2.6, index = [1, 1, 1])
    # cl = CutSurface(atoms, cutoffs = 1.2)
    atoms = cl.build()
    atoms = cl.rotate()
    print(atoms)
    view(atoms)

