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


class ConnectivityList:
    """
    Connectivity list object.

    threshold: float
        If the distance of two points smaller than the threshold, they will be counted as
        connected. The (Cartesian) distance is defined as: d(p, q) = sqrt((px − qx)**2 + (py − qy)**2).
    

    Example:
      cl = ConnectivityList(threshold)
      cl.build(points)
      output = cl.get_connectivity(pairs)
    """

    def __init__(self, cutoffs, search_list):
        self.cutoffs = cutoffs
        self.search_list = search_list
    def get_connectivity(self, inds):
        """
        pairs: list
            Pairs of points index

        Return:

        output: list
             A list of booleans. where the i−th element in the output list is true if
             the ind of points at the i−th position of the input pairs list are connected, 
             false otherwise.         
        """
        output = []
        groups = []
        for ind in inds:
            # print(ind)
            groups = groups + (self.groups[self.indexToGroup[ind]])
        return groups

    def update(self, iatom, groups):
        """
        Update the groups when adding a new point.

        Parameters:

        iatom: int
            The index of the point in the points array
        groups: list
            List contains all the groups. Points inside a group are connected.

        Returns:

        newgroup: list
            The updated groups after adding the new point.
        newboundarys: list
            The updated boundarys after adding the new point.
        """
        ngroup = len(groups)
        tempgroup = [iatom]
        ingroup = [True]*ngroup
        # Loop groups
        for i in range(ngroup):
            # Loop boundary point in the group i 
            for j in groups[i]:
                # Group i and iatom are connected, merge them to tempgroup
                dis = self.atoms[iatom].position - self.atoms[j].position
                cf = (self.cutoffs[iatom] + self.cutoffs[j])
                if abs(dis[0]) > cf: continue
                if abs(dis[1]) > cf: continue
                if abs(dis[2]) > cf: continue
                dis = dis[0]*dis[0] + dis[1]*dis[1] + dis[2]*dis[2]
                cf = cf*cf
                # print(dis, cf)
                if dis < cf:
                    ingroup[i] = False
                    tempgroup += groups[i]
                    break
        # Delete old group which connected with iatom
        newgroup = [groups[i] for i in range(ngroup) if ingroup[i]]
        # Add iatom group
        newgroup.append(tempgroup)
        return newgroup

    def build(self, atoms):
        """Build the the connectivity list and groups.
        Parameters:

        points: list
            A list of point which contains the x and y coordinates.
        """
        self.atoms = atoms
        self.natom = len(atoms)
        print('number of atoms: ', self.natom)
        self.cutoffs = self.cutoffs * covalent_radii[self.atoms.numbers]
        # Start with 1 group
        groups = [[0]]
        links = []
        for i in range(1, self.natom):
            groups = self.update(i, groups)
            # print(i, groups)
        self.groups = groups
        # indexToGroup
        ngroup = len(groups)
        indexToGroup = [0]*self.natom
        # print(self.natom)
        for i in range(ngroup):
            for j in groups[i]:
                indexToGroup[j] = i
        self.indexToGroup = indexToGroup


if __name__ == "__main__":
    from ase.io import read, write
    from ase.build import molecule
    from ase.visualize import view
    from ase.data import covalent_radii
    from blase.tools import write_blender, get_polyhedra_kinds, get_bondpairs, search_pbc
    atoms = read('../examples/perovskite.cif')
    # atoms = search_pbc(atoms, remove_dict = {'I': ['Pb']})
    atoms = molecule('CH3CH2NH2')
    atoms.center(vacuum=3.0)
    # print(atoms.cell)
    atoms = atoms*[1, 1, 2]
    view(atoms)
    cl = ConnectivityList(cutoffs = 1.2, search_list = ['C', 'N'])
    cl.build(atoms)
    print(cl.groups)
    inds = cl.get_connectivity([1, 10])
    print(inds)

