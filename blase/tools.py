import numpy as np
import pickle
import os
import json
from ase import Atoms, Atom
from ase.data import covalent_radii, atomic_numbers, chemical_symbols
from ase.visualize import view
import pprint
import time
import copy

def get_bondpairs_2(atoms, add_bonds = {}, remove_bonds = {}):
    """
    """
    from blase.default_data import default_bonds
    for ele1, ele2 in add_bonds:
        default_bonds[ele1] += ele2
    for ele1, ele2 in remove_bonds:
        if ele2 in default_bonds[ele1]:
            default_bonds[ele1].remove(ele2)

    symbols = set(atoms.get_chemical_symbols())
    bondparis = {}
    for symbol in symbols:
        bondparis[symbol] = []
        for symbol in symbols:
            if symbol in default_bonds[symbol]:
                bondparis[symbol].append(symbol)
    return bondpairs

def get_bondpairs(atoms, cutoff=1.0, add_bonds = {}, remove_bonds = {}):
    """
    The default bonds are stored in 'default_bonds'
    Get all pairs of bonding atoms
    remove_bonds
    """
    from ase.data import covalent_radii
    from ase.neighborlist import NeighborList, NewPrimitiveNeighborList
    from blase.default_data import default_bonds

    for ele1, ele2 in add_bonds:
        default_bonds[ele1] += ele2
    tstart = time.time()
    cutoffs = cutoff * covalent_radii[atoms.numbers]
    nl = NeighborList(cutoffs=cutoffs, self_interaction=False, bothways=True, primitive=NewPrimitiveNeighborList)
    nl.update(atoms)
    # bondpairs = []
    bondpairs = {}
    natoms = len(atoms)
    for a in range(natoms):
        bondpairs[a] = []
        indices, offsets = nl.get_neighbors(a)
        # print(a, indices)
        ele1 = atoms[a].symbol
        for a2, offset in zip(indices, offsets):
            ele2 = atoms[a2].symbol
            flag = False
            if ele2 in default_bonds[ele1] or ele1 in default_bonds[ele2]:
                flag = True
            for key, kinds in remove_bonds.items():
                for kind in kinds:
                    if atoms[a].symbol == key and kind == '*' \
                       or atoms[a].symbol == key and atoms[a2].symbol == kind \
                       or atoms[a].symbol == kind and atoms[a2].symbol == key:
                        flag = False
            if flag:
                bondpairs[a].append([a2, offset])
    print('get_bondpairs: {0:10.2f} s'.format(time.time() - tstart))
    return bondpairs



def default_atom_kind(element, positions, color = 'jmol'):
    """
    """
    from ase.data import chemical_symbols
    from ase.data.colors import jmol_colors, cpk_colors
    atom_kind = {}
    number = chemical_symbols.index(element)
    if color.upper() == 'JMOL':
        color = jmol_colors[number]
    elif color.upper() == 'CPK':
        color = jmol_colors[number]
    radius = covalent_radii[number]
    atom_kind['element'] = element
    atom_kind['color'] = color
    atom_kind['transmit'] = 1.0
    atom_kind['radius'] = radius
    atom_kind['positions'] = positions
    atom_kind['balltype'] = None
    return atom_kind
def default_bond_kind(element, color = 'jmol'):
    """
    """
    from ase.data import chemical_symbols
    from ase.data.colors import jmol_colors, cpk_colors
    bond_kind = {}
    number = chemical_symbols.index(element)
    if color.upper() == 'JMOL':
        color = jmol_colors[number]
    elif color.upper() == 'CPK':
        color = jmol_colors[number]
    bond_kind['element'] = element
    bond_kind['color'] = color
    bond_kind['transmit'] = 1.0
    bond_kind['lengths'] = []
    bond_kind['centers'] = []
    bond_kind['normals'] = []
    bond_kind['verts'] = []
    bond_kind['faces'] = []
    return bond_kind
def get_atom_kinds(atoms, scale = 1.0, props = {}):
    tstart = time.time()
    if 'kinds' not in atoms.info:
        atoms.info['kinds'] = atoms.get_chemical_symbols()
    kinds = list(set(atoms.info['kinds']))
    if isinstance(scale, float):
        newscale = {kind:[scale, scale, scale] for kind in kinds}
    elif isinstance(scale, dict):
        newscale = scale
    # print(kinds)
    atom_kinds = {}
    for kind in kinds:
        atom_kinds[kind] = {}
        element = kind.split('_')[0]
        inds = [atom.index for atom in atoms if atoms.info['kinds'][atom.index]==kind]
        atom_kind = default_atom_kind(element, atoms.positions[[inds]])
        atom_kinds[kind] = atom_kind
        atom_kind['radius'] = atom_kind['radius']
        atom_kind['scale'] = newscale[kind]
        
        if props:
            if kind in props.keys():
                for prop, value in props[kind].items():
                    atom_kinds[kind][prop] = value
    print('get_atom_kinds: {0:10.2f} s'.format(time.time() - tstart))
    return atom_kinds
def get_bond_kinds(atoms, bondlist):
    '''
    Build faces for instancing bonds.
    The radius of bonds is determined by nbins.
    mesh.from_pydata(vertices, [], faces)
    '''
    # view(atoms)
    
    tstart = time.time()
    # bond_kinds = copy.deepcopy(atom_kinds)
    bond_kinds = {}
    if 'kinds' not in atoms.info:
        atoms.info['kinds'] = atoms.get_chemical_symbols()
    for ind1, pairs in bondlist.items():
        kind1 = atoms.info['kinds'][ind1]
        element = kind1.split('_')[0]
        bond_kind = default_bond_kind(element)
        if kind1 not in bond_kinds:
            bond_kinds[kind1] = bond_kind
        # print(ind1, kind, pairs)
        for bond in pairs:
            ind2, offset = bond
            kind2 = atoms.info['kinds'][ind2]
            R = np.dot(offset, atoms.cell)
            pos = [atoms.positions[ind1],
                   atoms.positions[ind2] + R]
            center0 = (pos[0] + pos[1])/2.0
            vec = pos[0] - pos[1]
            length = np.linalg.norm(vec)
            nvec = vec/length
            nvec = nvec + 1e-8
            # verts, faces
            v1 = nvec + np.array([1.2323, 0.493749, 0.5604937284])
            v11 = v1 - np.dot(v1, nvec)*nvec
            v11 = v11/np.linalg.norm(v11)/2.828427
            v22 = np.cross(nvec, v11)*length*length
            #
            center = (center0 + pos[0])/2.0
            bond_kinds[kind1]['centers'].append(center)
            bond_kinds[kind1]['lengths'].append(length/4.0)
            bond_kinds[kind1]['normals'].append(nvec)
            nvert = len(bond_kinds[kind1]['verts'])
            bond_kinds[kind1]['verts'].append(center + v11)
            bond_kinds[kind1]['verts'].append(center - v11)
            bond_kinds[kind1]['verts'].append(center + v22)
            bond_kinds[kind1]['verts'].append(center - v22)
            bond_kinds[kind1]['faces'].append([nvert + 0, nvert + 2, nvert + 1, nvert + 3])
    # pprint.pprint(bond_kinds)
    print('get_bond_kinds: {0:10.2f} s'.format(time.time() - tstart))
    return bond_kinds

def get_polyhedra_kinds(atoms, atom_kinds, bondlist = {}, transmit = 0.8, polyhedra_dict = {}):
    """
    Two modes:
    (1) Search atoms bonded to kind
    polyhedra_dict: {'kind': ligands}
    """
    from scipy.spatial import ConvexHull
    from ase.data import covalent_radii
    from ase.neighborlist import NeighborList
    tstart = time.time()
    polyhedra_kinds = {}
    # loop center atoms
    # for ind1, pairs in bondlist.items():
        # kind = atoms.info['kinds'][ind1]
    for kind, ligand in polyhedra_dict.items():
        # print(kind, ligand)
        if kind not in polyhedra_kinds.keys():
            element = kind.split('_')[0]
            vertices = []
            edges = []
            faces = []
            polyhedra_kinds[kind] = {'vertices': vertices, 'edges': edges, 'faces': faces}
            lengths = []
            centers = []
            normals = []
            polyhedra_kinds[kind]['edge_cylinder'] = {'lengths': lengths, 'centers': centers, 'normals': normals}
            # number = chemical_symbols.index(kind)
            # color = jmol_colors[number]
            # polyhedra_kinds[kind]['number'] = number
            polyhedra_kinds[kind]['color'] = atom_kinds[kind]['color']
            polyhedra_kinds[kind]['transmit'] = transmit
            polyhedra_kinds[kind]['edge_cylinder']['color'] = (1.0, 1.0, 1.0)
            polyhedra_kinds[kind]['edge_cylinder']['transmit'] = 1.0
        inds = [atom.index for atom in atoms if atom.symbol == kind]
        for ind in inds:
            vertice = []
            for bond in bondlist[ind]:
            # indices, offsets = nl.get_neighbors(ind)
            # for a2, offset in zip(indices, offsets):
                a2, offset = bond
                if atoms[a2].symbol in ligand:
                    temp_pos = atoms[a2].position + np.dot(offset, atoms.cell)
                    vertice.append(temp_pos)
            nverts = len(vertice)
            # print(ind, indices, nverts)
            if nverts >3:
                # print(ind, vertice)
                # search convex polyhedra
                hull = ConvexHull(vertice)
                face = hull.simplices
                #
                # print(ind)
                nverts = len(polyhedra_kinds[kind]['vertices'])
                face = face + nverts
                edge = []
                for f in face:
                    edge.append([f[0], f[1]])
                    edge.append([f[0], f[2]])
                    edge.append([f[1], f[2]])
                polyhedra_kinds[kind]['vertices'] = polyhedra_kinds[kind]['vertices'] + list(vertice)
                polyhedra_kinds[kind]['edges'] = polyhedra_kinds[kind]['edges'] + list(edge)
                polyhedra_kinds[kind]['faces'] = polyhedra_kinds[kind]['faces'] + list(face)
                #
                # print('edge: ', edge)
                for e in edge:
                    # print(e)
                    center = (polyhedra_kinds[kind]['vertices'][e[0]] + polyhedra_kinds[kind]['vertices'][e[1]])/2.0
                    vec = polyhedra_kinds[kind]['vertices'][e[0]] - polyhedra_kinds[kind]['vertices'][e[1]]
                    length = np.linalg.norm(vec)
                    nvec = vec/length
                    # print(center, nvec, length)
                    polyhedra_kinds[kind]['edge_cylinder']['lengths'].append(length/2.0)
                    polyhedra_kinds[kind]['edge_cylinder']['centers'].append(center)
                    polyhedra_kinds[kind]['edge_cylinder']['normals'].append(nvec)
    print('get_polyhedra_kinds: {0:10.2f} s'.format(time.time() - tstart))
    return polyhedra_kinds


def euler_from_vector(normal, s = 'zxy'):
    from scipy.spatial.transform import Rotation as R
    normal = normal/np.linalg.norm(normal)
    vec = np.cross([0.0000014159, 0.000001951, 1], normal)
    vec = vec/np.linalg.norm(vec)
    # print(vec)
    # ang = np.arcsin(np.linalg.norm(vec))
    ang = np.arccos(normal[2])
    vec = -1*ang*vec
    # print(vec)
    r = R.from_rotvec(vec)
    euler = r.as_euler()
    return euler


def getEquidistantPoints(p1, p2, n):
    return zip(np.linspace(p1[0], p2[0], n+1), np.linspace(p1[1], p2[1], n+1), np.linspace(p1[2], p2[2], n+1))

def search_pbc(atoms, cutoff = [0.01, 0.01, 0.01]):
    """
    cutoffs: float or list
    """
    if isinstance(cutoff, float):
        cutoff = [cutoff]*3
    cutoff = 0.5 - np.array(cutoff)
    print('Search pbc: ')
    natoms = len(atoms)
    positions = atoms.get_scaled_positions()
    bdatoms = Atoms()
    index = {}
    na = 0
    for i in range(natoms):
        index[i] = []
        dv = 0.5 - positions[i]
        flag = abs(dv) > cutoff
        # print(dv, cutoff,  flag)
        flag = flag*1 + 1
        for l in range(flag[0]):
            for m in range(flag[1]):
                for n in range(flag[2]):
                    if l == 0 and m == 0 and n == 0: continue
                    temp_pos = positions[i] + np.sign(dv)*[l, m, n]
                    temp_pos = temp_pos.dot(atoms.cell)
                    bdatoms = bdatoms + Atom(atoms[i].symbol, temp_pos)
                    index[i].append(na)
                    na += 1

    return bdatoms, index

def get_cell_vertices(cell):
    cell_vertices = np.empty((2, 2, 2, 3))
    for c1 in range(2):
        for c2 in range(2):
            for c3 in range(2):
                cell_vertices[c1, c2, c3] = np.dot([c1, c2, c3],
                                                    cell)
    cell_vertices.shape = (8, 3)
    return cell_vertices
def get_bbox(bbox, atoms, show_unit_cell = True):
    """
    """
    if bbox is None:
        bbox = np.zeros([3, 2])
        positions = atoms.positions
        print(positions)
        cell_vertices = get_cell_vertices(atoms.cell)
        R = 4.0
        for i in range(3):
            P1 = (positions[:, i] - R).min(0)
            P2 = (positions[:, i] + R).max(0)
            if show_unit_cell and cell_vertices is not None:
                C1 = (cell_vertices[:, i]).min(0)
                C2 = (cell_vertices[:, i]).max(0)
                P1 = min(P1, C1)
                P2 = max(P2, C2)
            bbox[i] = [P1, P2]
        bbox = bbox
        w = bbox[0][1] - bbox[0][0]
        h = bbox[1][1] - bbox[1][0]
    else:
        w = (bbox[0][1] - bbox[0][0])
        h = (bbox[1][1] - bbox[1][0])
    return bbox, w, h


if __name__ == "__main__":
    from ase.io import read, write
    from ase.build import molecule, bulk
    from ase.visualize import view
    from ase.data import covalent_radii
    # test pbc
    atoms = bulk('Pt', cubic = True)
    # atoms = read('datas/ceo2.cif')
    atoms[0].z += 0.1
    bd_atoms, index = search_pbc(atoms, cutoff = [0.01, 0.01, 0.00])
    view(atoms + bd_atoms)
    print(bd_atoms.positions)
    print(index)
    # view(atoms)
    