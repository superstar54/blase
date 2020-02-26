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

    def __init__(self, atoms, cutoffs, search_dict = None, remove_dict = None):
        self.atoms = atoms
        self.natoms = len(atoms)
        self.cell = atoms.cell
        self.cutoffs = cutoffs
        self.remove_dict = remove_dict
        self.search_dict = search_dict
        self.add_list = {}
        # print(images)
        self.flags = {}
        self.origin = {}
        for i in range(self.natoms):
            self.flags[i] = ['000']
            self.origin[i] = [i, [0, 0, 0]]
        newatoms = Atoms()
        for x in [0, -1, 1]:
            for y in [0, -1, 1]:
                for z in [0, -1, 1]:
                    temp = self.atoms.copy()
                    t = x*self.cell[0] + y*self.cell[1] + z*self.cell[2]
                    # print(x, y, z, t)
                    temp.translate(t)
                    newatoms = newatoms + temp
        #
    def build(self, atoms = None, search_list = None):
        """
        """
        self.search_pbc()
        # view(atoms)
        # print(atoms)
        self.search_bond()
        # print(atoms)
        return self.atoms


    
    def search_molecule(self, atoms, cutoff=1.2, search_list = None):
        # atoms = atoms*[3, 3, 3]
        tstart = time.time()
        natoms = len(atoms)
        # print(natoms)
        cell = atoms.cell
        # view(atoms)
        if not search_list:
            formula = atoms.symbols.formula
            eles = formula.count()
            search_list = list(eles.keys())
        tstart = time.time()
        cl = ConnectivityList(cutoffs = 1.2, search_list = search_list)
        cl.build(atoms, search_list)
        print(cl.groups)
        print('time: {0:1.2f} s'.format(time.time() - tstart))
        mask = list(range(natoms))
        tstart = time.time()
        for i in range(natoms):
            # print(i)
            if atoms[i].symbol in search_list:
                for group in cl.groups:
                    if i in group:
                        mask = mask + group
        print('time: {0:1.2f} s'.format(time.time() - tstart))
        mask = list(set(mask) - set(range(natoms)))
        # print(mask)
        # view(atoms[mask])
        atoms = atoms + atoms[mask]
        return atoms
    
    def search_pbc(self, atoms = None, cutoff=1e-6):
        """
        Two modes:
        (1) Search atoms bonded to kind
        search_dict: {'kind': ligands}
        """
        from ase import Atom
        tstart = time.time()
        # loop center atoms
        if not atoms:
            atoms = self.atoms
        for i in range(3):
            natoms = len(atoms)
            positions = atoms.get_scaled_positions()
            for j in range(natoms):
                flag = 0.5 - positions[j, i]
                flag1 = abs(positions[j, i] % 1)
                if abs(flag1) <  cutoff:
                    flag = int(flag/abs(flag))
                    temp_pos = positions[j].dot(atoms.cell) + flag*atoms.cell[i]*0.99999
                    atoms = atoms + Atom(atoms[j].symbol, temp_pos)
                    ind = len(atoms) - 1
                    self.origin[ind] = [self.origin[j][0], self.origin[j][1].copy()]
                    self.origin[ind][1][i] += flag
                    if self.origin[ind][1] not in self.flags[self.origin[ind][0]]:
                        key = '{0}{1}{2}'.format(self.origin[ind][1][0], self.origin[ind][1][1], self.origin[ind][1][2])
                        self.flags[self.origin[ind][0]].append(key)
        print('search_pbc time: {0:1.2f} s'.format(time.time() - tstart))
        # print(self.origin)
        self.atoms = atoms
        return atoms
    def search_bond(self, atoms = None, cutoff=1.0):
        # search bonds
        # view(atoms)
        # search bonds out of unit cell
        tstart = time.time()
        if not atoms:
            atoms = self.atoms
        if not self.search_dict:
            self.search_dict = {}
            formula = atoms.symbols.formula
            eles = formula.count()
            keys = list(eles.keys())
            for ele in eles:
                self.search_dict[ele] = keys.copy()
                self.search_dict[ele].remove(ele)
        if self.remove_dict:
            for kind, eles in self.remove_dict.items():
                for ele in eles:
                    self.search_dict[kind].remove(ele)
        #
        newatoms = atoms.copy()
        cutoffs = cutoff * covalent_radii[atoms.numbers]
        nl = NeighborList(cutoffs=cutoffs, self_interaction=False, bothways=True)
        nl.update(atoms)
        for kind, ligand in self.search_dict.items():
            inds = [atom.index for atom in atoms if atom.symbol == kind]
            for ind in inds:
                indices, offsets = nl.get_neighbors(ind)
                for a2, offset in zip(indices, offsets):
                    flag = abs(offset[0]) + abs(offset[1]) + abs(offset[2])
                    if atoms[a2].symbol in ligand and flag > 0:
                        offset0 = [self.origin[a2][1][i] + offset[i] for i in range(3)]
                        key = '{0}{1}{2}'.format(offset0[0], offset0[1], offset0[2])
                        if key not in self.flags[self.origin[a2][0]]:
                            self.flags[self.origin[a2][0]].append(key)
                            temp_pos = atoms[a2].position + np.dot(offset, atoms.cell)
                            newatoms = newatoms + Atom(atoms[a2].symbol, temp_pos)
                            ind = len(newatoms) - 1
                            self.origin[ind] = [self.origin[a2][0], [self.origin[a2][1][i] + offset[i] for i in range(3)]]
        #
        print('search_bond time: {0:1.2f} s'.format(time.time() - tstart))
        # newatoms.pbc = [False, False, False]
        self.atoms = newatoms

        return newatoms


if __name__ == "__main__":
    from ase.io import read, write
    from ase.build import molecule, bulk
    from ase.visualize import view
    from ase.data import covalent_radii
    # atoms = read('../examples/tio2.cif')
    atoms = read('../examples/perovskite.cif')
    # view(atoms)
    cl = ConnectivityList(atoms, cutoffs = 1.2, search_dict = {'Ti':['O']})
    atoms = cl.build()
    view(atoms)
    # inds = cl.get_connectivity([1, 10])

