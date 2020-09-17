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

    def __init__(self, atoms, cutoffs, pbc = [], bonds_dict = {}, molecule_list = []):
        self.atoms = atoms
        self.natoms = len(atoms)
        self.cell = atoms.cell
        self.cutoffs = cutoffs
        self.pbc = pbc
        self.molecule_list = []
        self.molecule_dict = {}
        self.flags = {}
        self.origin = {}
        self.add_total = 0
        for i in range(self.natoms):
            self.flags[i] = ['000']
            self.origin[i] = [i, [0, 0, 0]]
        #
        self.bonds_dict = {}
        formula = self.atoms.symbols.formula
        eles = formula.count()
        keys = list(eles.keys())
        for ele in eles:
            ligands = keys.copy()
            ligands.remove(ele)
            self.bonds_dict[ele] = ligands
            self.molecule_dict[ele] = []
        for key, value in bonds_dict.items():
            # print(key, value)
            if value[1] == 0:
                self.bonds_dict[key] = value
            if value[1] == -1:
                for ele in value[0]:
                    self.bonds_dict[key].remove(ele)
        for eles in molecule_list:
            # print(eles)
            # print(self.molecule_dict[eles[0]])
            self.molecule_dict[eles[0]] = self.molecule_dict[eles[0]] + [eles[1]]
            self.molecule_dict[eles[1]] = self.molecule_dict[eles[1]] + [eles[0]]
            # self.molecule_dict[value[0]] += self.molecule_dict + [ele] + value[1:]
        print(self.bonds_dict)
        print(self.molecule_dict)
        #
    def build(self, atoms = None, pbc = None, bonds_dict = None, molecule_dict = None):
        """
        """
        self.add = 0
        tstart = time.time()
        self.search_pbc()
        # view(self.atoms)
        print('Start searching: ')
        self.search_bond(bonds_dict = self.bonds_dict)
        self.search_molecule(molecule_dict = self.molecule_dict)
        self.search_bond(bonds_dict = self.bonds_dict)
        print('Search ConnectivityList time: {0:1.2f} s, add {1} atoms'.format(time.time() - tstart, self.add_total))
        self.atoms.pbc = [False, False, False]
        print(self.atoms)
        return self.atoms
    def search_molecule(self, molecule_dict):
        # print('Search molecule: ')
        # atoms = atoms*[3, 3, 3]
        flag = True
        i = 0
        while flag:
            self.search_bond(bonds_dict = molecule_dict)
            i += 1
            if self.add == 0: flag = False
        # view(atoms[mask])
    def search_pbc(self, atoms = None, cutoff=1e-2):
        """
        Two modes:
        (1) Search atoms bonded to kind
        bonds_dict: {'kind': ligands}
        """
        # print('Search pbc: ')
        tstart = time.time()
        # loop center atoms
        if not atoms:
            atoms = self.atoms
        for i in range(3):
            natoms = len(atoms)
            positions = atoms.get_scaled_positions()
            for j in range(natoms):
                if atoms[j].symbol not in self.pbc: continue
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
        # print('search_pbc time: {0:1.2f} s'.format(time.time() - tstart))
        # print(self.origin)
        self.atoms = atoms
        return atoms
    def search_bond(self, atoms = None, cutoff=1.0, bonds_dict = None):
        # search bonds
        # view(atoms)
        # search bonds out of unit cell
        # print('Search bonds: ')
        # print(bonds_dict)
        tstart = time.time()
        if not atoms:
            atoms = self.atoms
        #
        self.add = 0
        newatoms = atoms.copy()
        cutoffs = cutoff * covalent_radii[atoms.numbers]
        nl = NeighborList(cutoffs=cutoffs, self_interaction=False, bothways=True)
        nl.update(atoms)
        for kind, ligand in bonds_dict.items():
            inds = [atom.index for atom in atoms if atom.symbol == kind]
            for ind in inds:
                indices, offsets = nl.get_neighbors(ind)
                # print(kind, ind, indices)
                for a2, offset in zip(indices, offsets):
                    flag = abs(offset[0]) + abs(offset[1]) + abs(offset[2])
                    if atoms[a2].symbol in ligand and flag > 0:
                        offset0 = [self.origin[a2][1][i] + offset[i] for i in range(3)]
                        key = '{0}{1}{2}'.format(offset0[0], offset0[1], offset0[2])
                        # print(a2, key)
                        if key not in self.flags[self.origin[a2][0]]:
                            self.flags[self.origin[a2][0]].append(key)
                            temp_pos = atoms[a2].position + np.dot(offset, atoms.cell)
                            newatoms = newatoms + Atom(atoms[a2].symbol, temp_pos)
                            ind = len(newatoms) - 1
                            self.origin[ind] = [self.origin[a2][0], [self.origin[a2][1][i] + offset[i] for i in range(3)]]
                            self.add += 1
            # print(kind, ligand, self.add)
        #
        # print('search_bond time: {0:1.2f} s'.format(time.time() - tstart))
        # newatoms.pbc = [False, False, False]
        # print('add atoms: ', self.add)
        self.atoms = newatoms
        self.add_total += self.add
        return newatoms
    def remove_pbc(self, atoms = None, cutoff=1.0, eps=1e-3):
        """
        Two modes:
        (1) Search atoms bonded to kind
        bonds_dict: {'kind': ligands}
        """
        # print('Search pbc: ')
        tstart = time.time()
        # loop center atoms
        if not atoms:
            atoms = self.atoms
        natoms = len(atoms)
        positions = atoms.get_scaled_positions()
        mask = []
        newpositions = []
        for i in range(natoms):
            if positions[i, 0] < 0 or \
               positions[i, 0] > 1 or \
               positions[i, 1] < 0 or \
               positions[i, 1] > 1 or \
               positions[i, 2] < 0 or \
               positions[i, 2] > 1:
                mask.append(i)
                continue


        del atoms[mask]
        #atoms.wrap()
        natoms = len(atoms)
        mask = []
        positions = atoms.get_scaled_positions()
        for i in range(natoms):
            if positions[i, 0] > 1 - eps: positions[i, 0] = 0
            if positions[i, 1] > 1 - eps: positions[i, 1] = 0
            if positions[i, 2] > 1 - eps: positions[i, 2] = 0
        for i in range(natoms - 1):
            for j in range(i + 1, natoms):
                pos = positions[i] - positions[j]
                dis = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]
                #dis = atoms.get_distance(i, j, mic=False)
                #print(i, j, dis)
                if dis < eps*eps:
                    print(i, j, dis)
                    if pos[0] < eps or pos[1] < eps or pos[2] < eps:
                        mask.append(j)
                    else:
                        mask.append(i)
                        break
        #
        #print(mask)
        del atoms[mask]
        self.atoms = atoms
        return atoms


if __name__ == "__main__":
    from ase.io import read, write
    from ase.build import molecule, bulk
    from ase.visualize import view
    from ase.data import covalent_radii
    # test pbc
    # atoms = bulk('Pt')
    # cl = ConnectivityList(atoms, cutoffs = 1.2, pbc = ['Pt'])
    # atoms = cl.build()
    # view(atoms)
    # atoms = read('../examples/datas/tio2.cif')
    # view(atoms)
    # cl = ConnectivityList(atoms, cutoffs = 1.2, pbc = ['Ti'], bonds_dict = {'O': [['Ti'], -1]})
    # atoms = read('../examples/datas/latino2-trans-p-001-plp-relax.in')
    # view(atoms)
    # cl = ConnectivityList(atoms, cutoffs = 1.2, pbc = ['La'], bonds_dict = {}, molecule_list = {})
    # atoms = read('../examples/perovskite.cif')
    # cl = ConnectivityList(atoms, cutoffs = 1.2,  bonds_dict = {'I': [['Pb'], -1]}, molecule_list = [['C', 'N']])
    # atoms = read('../examples/datas/anthraquinone.cif')
    # cl = ConnectivityList(atoms, cutoffs = 1.2, molecule_list = [['C', 'C'], ['C', 'O']])
    # print(atoms)
    # cl = ConnectivityList(atoms, cutoffs = 1.2)
    # atoms = cl.build()
    # print(atoms)
    # view(atoms)
    # atoms = cl.remove_pbc(atoms)
    # print(atoms)
    # view(atoms)

