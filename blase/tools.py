import numpy as np
import pickle
import os
import json
from ase import Atoms, Atom
from ase.data import covalent_radii, atomic_numbers, chemical_symbols
from ase.data.colors import jmol_colors
from ase.visualize import view
import pprint
import time
import copy

def get_bondpairs(atoms, cutoff=1.0, rmbonds = {}):
    """
    Get all pairs of bonding atoms
    rmbonds
    """
    from ase.data import covalent_radii
    from ase.neighborlist import NeighborList, NewPrimitiveNeighborList
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
        for a2, offset in zip(indices, offsets):
            flag = True
            for key, kinds in rmbonds.items():
                for kind in kinds:
                    if atoms[a].symbol == key and kind == '*' \
                       or atoms[a].symbol == key and atoms[a2].symbol == kind \
                       or atoms[a].symbol == kind and atoms[a2].symbol == key:
                        flag = False
            if flag:
                bondpairs[a].append([a2, offset])
    print('get_bondpairs: {0:10.2f} s'.format(time.time() - tstart))
    return bondpairs


def write_blender(atoms, display = False, queue = None, **kwargs):
    with open('blase.inp', 'wb') as f:
        pickle.dump([atoms, kwargs], f)
    #
    blender_cmd = 'blender'
    if 'BLENDER_COMMAND' in os.environ.keys():
        blender_cmd = os.environ['BLENDER_COMMAND']
    blase_path = os.environ['BLASE_PATH']
    blase_cmd = blase_path + '/bin/run-blase.py'
    if display:
        cmd = blender_cmd + ' -P ' + blase_cmd
    elif queue == 'SLURM':
        cmd = 'srun -n $SLURM_NTASKS ' +  blender_cmd + ' -b ' + ' -P ' + blase_cmd
    else:
        cmd = blender_cmd + ' -b ' + ' -P ' + blase_cmd
    print(cmd)
    errcode = os.system(cmd)
    # if errcode != 0:
    #     raise OSError('Command ' + cmd +
    #                   ' failed with error code %d' % errcode)

# def get_atom_kinds(atoms, props):
    # return kinds
def get_atom_kinds(atoms, props = {}):
    # symbols = atoms.symbols
    # formula = atoms.symbols.formula
    # atom_kinds = formula.count()
    tstart = time.time()
    if hasattr(atoms, 'kinds'):
        kinds = list(set(atoms.kinds))
    else:
        atoms.kinds = atoms.get_chemical_symbols()
        kinds = list(set(atoms.kinds))
    # print(kinds)
    atom_kinds = {}
    for kind in kinds:
        atom_kinds[kind] = {}
        element = kind.split('_')[0]
        print(kind, element)
        number = chemical_symbols.index(element)
        inds = [atom.index for atom in atoms if atoms.kinds[atom.index]==kind]
        color = jmol_colors[number]
        radius = covalent_radii[number]
        atom_kinds[kind]['element'] = element
        atom_kinds[kind]['positions'] = atoms[inds].positions
        atom_kinds[kind]['number'] = number
        atom_kinds[kind]['color'] = color
        atom_kinds[kind]['transmit'] = 1.0
        atom_kinds[kind]['radius'] = radius
        atom_kinds[kind]['balltype'] = None
        # bond
        atom_kinds[kind]['lengths'] = []
        atom_kinds[kind]['centers'] = []
        atom_kinds[kind]['normals'] = []
        atom_kinds[kind]['verts'] = []
        atom_kinds[kind]['faces'] = []
        if props:
            if kind in props.keys():
                for prop, value in props[kind].items():
                    atom_kinds[kind][prop] = value
    print('get_atom_kinds: {0:10.2f} s'.format(time.time() - tstart))
    return atom_kinds
def get_bond_kinds(atoms, atom_kinds, bondlist):
    '''
    Build faces for instancing bonds.
    The radius of bonds is determined by nbins.
    mesh.from_pydata(vertices, [], faces)
    '''
    # view(atoms)
    
    tstart = time.time()
    # bond_kinds = copy.deepcopy(atom_kinds)
    bond_kinds = atom_kinds.copy()
    for ind1, pairs in bondlist.items():
        kind1 = atoms.kinds[ind1]
        # print(ind1, kind, pairs)
        for bond in pairs:
            ind2, offset = bond
            kind2 = atoms.kinds[ind2]
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
            kinds = [kind1, kind2]
            for i in range(1):
                kind = kinds[i]
                center = (center0 + pos[i])/2.0
                bond_kinds[kind]['centers'].append(center)
                bond_kinds[kind]['lengths'].append(length/4.0)
                bond_kinds[kind]['normals'].append(nvec)
                nvert = len(bond_kinds[kind]['verts'])
                bond_kinds[kind]['verts'].append(center + v11)
                bond_kinds[kind]['verts'].append(center - v11)
                bond_kinds[kind]['verts'].append(center + v22)
                bond_kinds[kind]['verts'].append(center - v22)
                bond_kinds[kind]['faces'].append([nvert + 0, nvert + 2, nvert + 1, nvert + 3])
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
        # kind = atoms.kinds[ind1]
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

