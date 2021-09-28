import numpy as np
from ase import Atoms, Atom, atom
from ase.data import covalent_radii, atomic_numbers, chemical_symbols
from ase.visualize import view
import time
from ase.neighborlist import neighbor_list
from pprint import pprint

from numpy.lib.utils import source

def get_bondlists(atoms, bondsetting):
    """
    The default bonds are stored in 'default_bonds'
    Get all pairs of bonding atoms
    remove_bonds
    """
    from neighborlist import neighbor_list
    if len(bondsetting) == 0: return {}
    cutoff_min = {}
    cutoff = {}
    for key, data in bondsetting.items():
        cutoff_min[key] = data[0]
        cutoff[key] = data[1]
    #
    tstart = time.time()
    nli, nlj, nlS = neighbor_list('ijS', atoms, cutoff, cutoff_min, self_interaction=False)
    bondlists = np.append(np.array([nli, nlj], dtype=int).T, np.array(nlS, dtype=int), axis = 1)
    # inds = list(set(nli))
    # for i in inds:
        # bondlists[i].append([j, offset[0], offset[1], offset[2]])
    # bondlists = {i: data  for i, data in bondlists.items() if len(data) > 0}
    print('get_bondlists 1: {0:10.2f} s'.format(time.time() - tstart))
    return bondlists

def search_skin(atoms, bondsetting, bondlists, skin = []):
    """
    The default bonds are stored in 'default_bonds'
    Get all pairs of bonding atoms
    remove_bonds
    """
    from neighborlist import neighbor_list
    tstart = time.time()
    atoms_skin = Atoms()
    atoms_skin.info['species'] = []
    if len(bondsetting) == 0: return atoms_skin
    if 'species' not in atoms.info:
        atoms.info['species'] = atoms.get_chemical_symbols()
    speciesarray = np.array(atoms.info['species'])
    specieslist = list(set(atoms.info['species']))
    ncore = len(atoms) - len(skin)
    ind1 = bondlists[:, 0] < ncore
    ind2 = bondlists[:, 1] < ncore
    ind3 = bondlists[:, 0] > ncore
    ind4 = bondlists[:, 1] > ncore
    bondlists0 = bondlists[ind1&ind2]
    bondlistsij = bondlists[ind1&ind4]
    bondlistsji = bondlists[ind2&ind3]
    for spi in specieslist:
        bondlists1 = bondlistsij[speciesarray[bondlistsij[:, 0]] == spi]
        for spj in specieslist:
            if (spi, spj) not in bondsetting: continue
            if bondsetting[(spi, spj)][3]:
                bondlists0 = np.append(bondlists0, bondlists1[speciesarray[bondlists1[:, 1]] == spj], axis = 0)
        bondlists1 = bondlistsji[speciesarray[bondlistsij[:, 0]] == spi]
        for spj in specieslist:
            if (spj, spi) not in bondsetting: continue
            if bondsetting[(spj, spi)][3]:
                bondlists0 = np.append(bondlists0, bondlists1[speciesarray[bondlists1[:, 1]] == spj], axis = 0)
    # print(bondlists0)
    ind_atom_skin = list(set(bondlists0[:, 1]) & set(skin))
    if len(ind_atom_skin) == 0:
        return atoms_skin
    atoms_skin = atoms[ind_atom_skin]
    atoms_skin.info['species'] = speciesarray[ind_atom_skin]
    print('search skin : {0:10.2f} s'.format(time.time() - tstart))
    return atoms_skin

def calc_bond_data(atoms, bondlists):
    """
    """
    from tools import get_bond_kind
    from ase.data import chemical_symbols, covalent_radii
    chemical_symbols = np.array(chemical_symbols)
    if 'species' not in atoms.info:
        atoms.info['species'] = atoms.get_chemical_symbols()
    speciesarray = np.array(atoms.info['species'])
    elementarray = np.array(atoms.get_chemical_symbols())
    radiusarray = np.array([covalent_radii[chemical_symbols == em][0] for em in elementarray])
    specieslist = list(set(atoms.info['species']))
    bond_kinds = {}
    for spi in specieslist:
        bondlists1 = bondlists[speciesarray[bondlists[:, 0]] == spi]
        emi = spi.split('_')[0]
        emj = elementarray[bondlists1[:, 1]]
        bond_kind = get_bond_kind(emi)
        if spi not in bond_kinds:
            bond_kinds[spi] = bond_kind
        offset = bondlists1[:, 2:5]
        R = np.dot(offset, atoms.cell)
        vec = atoms.positions[bondlists1[:, 0]] - (atoms.positions[bondlists1[:, 1]] + R)
        length = np.linalg.norm(vec, axis = 1)
        nvec = vec/length[:, None]
        radius1 = covalent_radii[chemical_symbols ==emi][0]
        radius2 = radiusarray[bondlists1[:, 1]]
        pos = [atoms.positions[bondlists1[:, 0]],
                atoms.positions[bondlists1[:, 1]] + R ]
        center0 = (pos[0] + pos[1])/2.0
        vec = pos[0] - pos[1]
        length = np.linalg.norm(vec, axis = 1)
        nvec = vec/length[:, None]
        nvec = nvec + 1e-8
        # verts, faces
        v1 = nvec + np.array([1.2323, 0.493749, 0.5604937284])
        tempv = np.einsum("ij, ij->i", v1, nvec)
        v11 = v1 - (nvec.T*tempv).T
        templ = np.linalg.norm(v11, axis = 1)
        v11 = v11/templ[:, None]/2.828427
        tempv = np.cross(nvec, v11)
        v22 = (tempv.T*(length*length)).T
        #
        center = (center0 + pos[0])/2.0
        bond_kinds[spi]['centers'] = center
        bond_kinds[spi]['lengths'] = length/4.0
        bond_kinds[spi]['normals'] = nvec
        nvert = np.arange(len(center))
        bond_kinds[spi]['verts'] = center + v11
        bond_kinds[spi]['verts'] = np.append(bond_kinds[spi]['verts'], center - v11, axis = 0)
        bond_kinds[spi]['verts'] = np.append(bond_kinds[spi]['verts'], center + v22, axis = 0)
        bond_kinds[spi]['verts'] = np.append(bond_kinds[spi]['verts'], center - v22, axis = 0)
        bond_kinds[spi]['faces'] = [nvert + 0, nvert + 2, nvert + 1, nvert + 3]
    # pprint(bond_kinds)
    return bond_kinds

def cylinder_mesh_from_instance(centers, normals, lengths, scale, source):
    from scipy.spatial.transform import Rotation as R
    
    tstart = time.time()
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
        ang = np.arccos(normal[2]*0.999999)
        vec = -1*ang*vec
        r = R.from_rotvec(vec)
        matrix = r.as_matrix()
        # print(vec, ang)
        vert1 = vert0.copy()
        vert1 = vert1*np.array([scale, scale, length])
        vert1 = vert1.dot(matrix)
        vert1 += center
        for vert in vert1:
            verts.append(vert)
        for face in face0:
            face = [x+ i*nvert for x in face]
            faces.append(face)
    print('cylinder_mesh_from_instance: {0:10.2f} s'.format( time.time() - tstart))
    return verts, faces

def cylinder_mesh_from_instance_vec(centers, normals, lengths, scale, source):
    from scipy.spatial.transform import Rotation as R
    tstart = time.time()
    verts = []
    faces = []
    vert0, face0 = source
    vert0 = np.array(vert0)
    nvert = len(vert0)
    scale = np.array([[scale, scale]]*len(centers))
    print(scale)
    print(lengths)
    scale = np.append(scale, lengths.reshape(-1, 1), axis = 1)
    # scale = np.tile(scale, (len(vert0), 1))
    # normals = normals/np.linalg.norm(normals)
    print(np.cross([0.0000014159, 0.000001951, 1], normals[0]))
    vec = np.cross([0.0000014159, 0.000001951, 1], normals)
    # print(center, normals, vec)
    vec = vec/np.linalg.norm(vec, axis = 1)[:, None]
    # print(vec)
    print(np.arccos(normals[0, 2]*0.999999))
    ang = np.arccos(normals[:, 2]*0.999999)
    print(-1*ang[0]*vec[0])
    vec = -1*(vec.T*ang).T
    print(R.from_rotvec(vec[0]).as_matrix())
    r = R.from_rotvec(vec)
    matrix = r.as_matrix()
    # print(vec, ang)
    vert1 = vert0.copy()
    vert1 = [vert0]*len(centers)
    vert1 = vert1*scale
    vert1 = vert1.dot(matrix)
    vert1 += centers
    verts.append(vert1)
    # faces.append(face)
    print('cylinder_mesh_from_instance: {0:10.2f} s'.format( time.time() - tstart))
    return verts, faces



if __name__ == "__main__":
    from ase.build import bulk, molecule
    from ase.atoms import Atoms
    from ase.io import read
    from ase.visualize import view
    from ase.neighborlist import neighbor_list
    from boundary import search_boundary
    # atoms = bulk('Pt', cubic = True)
    atoms = read('docs/source/_static/datas/tio2.cif')
    # atoms = read('docs/source/_static/datas/mof-5.cif')
    # atoms = atoms*[4, 4, 4]
    # positions1, offsets1, positions2, offsets2 = search_boundary(atoms.positions, atoms.cell, boundary=[[-0.6, 1.6], [-0.6, 1.6], [-0.6, 1.6]])
    bondsetting = {('Ti', 'O'): [0, 2.5, False, False], ('O', 'O'): [0, 1.5, False, False]}
#     bondsetting = {('H' ,   'O' ):    [0.500,   1.164,      False,        False      ],
# ('O' ,   'H' ):    [0.500,   1.164,      False,        False      ],
# ('O' ,   'O' ):    [0.500,   1.584,      False,        False      ],
# ('Zn',    'H'):     [0.500,   1.836,      False,        False      ],
# ('Zn',    'O'):     [0.000,   2.500,      True ,        False      ],
# ('C' ,   'H' ):    [0.000,   1.400,      False,        False      ],
# ('C' ,   'O' ):    [0.500,   1.704,      False,        False      ],
# ('C' ,   'C' ):    [0.500,   1.824,      False,        False],}
    atoms = molecule('H2O')
    bondsetting = {('H', 'O'): [0, 1.5, False, False]}
    bondlists = get_bondlists(atoms, bondsetting)
    datas = calc_bond_data(atoms, bondlists)
    source = ([[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]], [[0, 1, 2, 3]])
    cylinder_mesh_from_instance(datas['O']['centers'], datas['O']['normals'], datas['O']['lengths'], 0.1, source)
    # print(bondlists)
    