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

def get_bondpairs(atoms, cutoff=1.2, add_bonds = {}, remove_bonds = {}):
    """
    The default bonds are stored in 'default_bonds'
    Get all pairs of bonding atoms
    remove_bonds
    """
    from ase.data import covalent_radii
    from ase.neighborlist import neighbor_list
    from blase.default_data import default_bonds
    from ase.neighborlist import natural_cutoffs

    for ele1, ele2 in add_bonds:
        default_bonds[ele1] += ele2
    tstart = time.time()
    if isinstance(cutoff, float):
        cutoff = cutoff * covalent_radii[atoms.numbers]
    # print(cutoff)
    nli, nlj, nlS = neighbor_list('ijS', atoms, cutoff=cutoff, self_interaction=False)
    bondpairs = {i: [] for i in set(nli)}
    if 'species' not in atoms.info:
        atoms.info['species'] = atoms.get_chemical_symbols()
    # print(atoms.info['species'])
    for i, j, offset in zip(nli, nlj, nlS):
        # i, j, offset = nl
        ele1 = atoms.info['species'][i]
        ele2 = atoms.info['species'][j]
        flag = False
        if ele2.split('_')[0] in default_bonds[ele1.split('_')[0]] or ele1.split('_')[0] in default_bonds[ele2.split('_')[0]]:
            flag = True
        for key, kinds in remove_bonds.items():
            for kind in kinds:
                # print(key, kind, atoms.info['species'][i], atoms.info['species'][j])
                if atoms.info['species'][i] == key and kind == '*' \
                    or atoms.info['species'][i] == key and atoms.info['species'][j] == kind \
                    or atoms.info['species'][i] == kind and atoms.info['species'][j] == key:
                    flag = False
                    # print(flag)
        if flag:
            bondpairs[i].append([j, offset])
    print('get_bondpairs: {0:10.2f} s'.format(time.time() - tstart))
    return bondpairs

def default_atom_kind(element, positions, color = 'JMOL'):
    """
    """
    from ase.data import chemical_symbols
    from ase.data.colors import jmol_colors, cpk_colors
    from blase.default_data import vesta_color
    atom_kind = {}
    number = chemical_symbols.index(element)
    if color.upper() == 'JMOL':
        color = jmol_colors[number]
    elif color.upper() == 'CPK':
        color = jmol_colors[number]
    elif color.upper() == 'VESTA':
        color = vesta_color[element]
    radius = covalent_radii[number]
    atom_kind['element'] = element
    atom_kind['color'] = color
    atom_kind['transmit'] = 1.0
    atom_kind['radius'] = radius
    atom_kind['scale'] = [1, 1, 1]
    atom_kind['positions'] = positions
    atom_kind['balltype'] = None
    return atom_kind
def default_bond_kind(element, color = 'JMOL'):
    """
    """
    from ase.data import chemical_symbols
    from ase.data.colors import jmol_colors, cpk_colors
    from blase.default_data import vesta_color
    # print('color: ', color)
    bond_kind = {}
    number = chemical_symbols.index(element)
    if color.upper() == 'JMOL':
        color = jmol_colors[number]
    elif color.upper() == 'CPK':
        color = jmol_colors[number]
    elif color.upper() == 'VESTA':
        color = vesta_color[element]
    bond_kind['element'] = element
    bond_kind['color'] = color
    bond_kind['transmit'] = 1.0
    bond_kind['lengths'] = []
    bond_kind['centers'] = []
    bond_kind['normals'] = []
    bond_kind['verts'] = []
    bond_kind['faces'] = []
    return bond_kind
def get_atom_kinds(atoms, scale = 1.0, props = {}, color = 'JMOL'):
    tstart = time.time()
    if 'species' not in atoms.info:
        atoms.info['species'] = atoms.get_chemical_symbols()
    kinds = list(set(atoms.info['species']))
    if isinstance(scale, float):
        newscale = {kind:[scale, scale, scale] for kind in kinds}
    elif isinstance(scale, dict):
        newscale = scale
    # print(kinds)
    atom_kinds = {}
    # print(atoms.info['species'])
    for kind in kinds:
        print(kind)
        atom_kinds[kind] = {}
        element = kind.split('_')[0]
        print(len(atoms), len(atoms.info['species']))
        inds = [atom.index for atom in atoms if atoms.info['species'][atom.index]==kind]
        atom_kind = default_atom_kind(element, atoms.positions[inds], color = color)
        atom_kinds[kind] = atom_kind
        atom_kind['radius'] = atom_kind['radius']
        atom_kind['scale'] = newscale[kind]
        
        if props:
            if kind in props.keys():
                for prop, value in props[kind].items():
                    atom_kinds[kind][prop] = value
    print('get_atom_kinds: {0:10.2f} s'.format(time.time() - tstart))
    return atom_kinds
def get_bond_kinds(atoms, atoms_kinds, bondlist, color = 'JMOL'):
    '''
    Build faces for instancing bonds.
    The radius of bonds is determined by nbins.
    mesh.from_pydata(vertices, [], faces)
    '''
    # view(atoms)
    
    tstart = time.time()
    # bond_kinds = copy.deepcopy(atom_kinds)
    bond_kinds = {}
    if 'species' not in atoms.info:
        atoms.info['species'] = atoms.get_chemical_symbols()
    for ind1, pairs in bondlist.items():
        kind1 = atoms.info['species'][ind1]
        element = kind1.split('_')[0]
        bond_kind = default_bond_kind(element, color = color)
        if kind1 not in bond_kinds:
            bond_kinds[kind1] = bond_kind
        # print(ind1, kind, pairs)
        for bond in pairs:
            ind2, offset = bond
            kind2 = atoms.info['species'][ind2]
            R = np.dot(offset, atoms.cell)
            vec = atoms.positions[ind1] - (atoms.positions[ind2] + R)
            length = np.linalg.norm(vec)
            nvec = vec/length
            pos = [atoms.positions[ind1] - nvec*atoms_kinds[kind1]['radius']*atoms_kinds[kind1]['scale'][0]*0.5,
                   atoms.positions[ind2] + R + nvec*atoms_kinds[kind2]['radius']*atoms_kinds[kind2]['scale'][0]*0.5]
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

def get_polyhedras(atoms, atom_kinds, bondlist = {}, polyhedra_dict = {}):
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
        # kind = atoms.info['species'][ind1]
    for kind, ligand in polyhedra_dict.items():
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
        # kind = atoms.info['species'][ind1]
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
    bdatoms = Atoms()
    bdatoms.info['species'] = []
    if 'species' not in atoms.info:
        atoms.info['species'] = atoms.get_chemical_symbols()
    index = {}
    if isinstance(cutoff, float):
        cutoff = [cutoff]*3
    cutoff = 0.5 - np.array(cutoff)
    print('Search pbc: ')
    natoms = len(atoms)
    positions = atoms.get_scaled_positions()
    
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
                    spe = bdatoms.info['species']
                    bdatoms = bdatoms + Atom(atoms[i].symbol, temp_pos)
                    spe.append(atoms.info['species'][i])
                    bdatoms.info['species'] = spe
                    index[i].append(na)
                    na += 1
    return bdatoms, index

def get_cell_vertices(cell):
    """
    """
    cell = np.array(cell)
    cell = cell.reshape(3, 3)
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
        if atoms.pbc.any():
            cell_vertices = get_cell_vertices(atoms.cell)
        else:
            cell_vertices = None
        # print('cell_vertices: ', cell_vertices)
        R = 1
        for i in range(3):
            P1 = (positions[:, i] - R).min(0)
            P2 = (positions[:, i] + R).max(0)
            # print('P1, P2: ', P1, P2)
            if show_unit_cell and cell_vertices is not None:
                C1 = (cell_vertices[:, i]).min(0)
                C2 = (cell_vertices[:, i]).max(0)
                P1 = min(P1, C1)
                P2 = max(P2, C2)
                print('C1, C2: ', C1, C2)
            bbox[i] = [P1, P2]
        bbox = bbox
    return bbox


if __name__ == "__main__":
    from ase.io import read, write
    from ase.build import molecule, bulk
    from ase.visualize import view
    from ase.data import covalent_radii
    # test bbox
    from ase.neighborlist import neighbor_list
    atoms = molecule('CH3OH')
    
    i = neighbor_list('ij', atoms, cutoff=1.0)
    print(i)
    view(atoms)
    