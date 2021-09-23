from operator import pos
import numpy as np
from ase import Atoms, Atom, atom
from ase.data import covalent_radii, atomic_numbers, chemical_symbols
from ase.visualize import view
import time



def get_bondpairs(atoms, bondsetting, skin = []):
    """
    The default bonds are stored in 'default_bonds'
    Get all pairs of bonding atoms
    remove_bonds
    """
    from ase.neighborlist import neighbor_list
    atoms_skin = Atoms()
    atoms_skin.info['species'] = []
    if len(bondsetting) == 0: return {}, atoms_skin
    cutoff_min = {}
    cutoff_max = {}
    for key, data in bondsetting.items():
        cutoff_min[key] = data[0]
        cutoff_max[key] = data[1]
    #
    tstart = time.time()
    nli_min, nlj_min, nlS_min = neighbor_list('ijS', atoms, cutoff=cutoff_min, self_interaction=False)
    nli_max, nlj_max, nlS_max = neighbor_list('ijS', atoms, cutoff=cutoff_max, self_interaction=False)
    bondpairs_min = {i: [] for i in set(nli_min)}
    bondpairs_max = {i: [] for i in set(nli_max)}
    print('get_bondpairs 1: {0:10.2f} s'.format(time.time() - tstart))
    tstart = time.time()
    if 'species' not in atoms.info:
        atoms.info['species'] = atoms.get_chemical_symbols()
    for i, j, offset in zip(nli_min, nlj_min, nlS_min):
        if i not in skin and j not in skin:
            bondpairs_min[i].append([j, offset[0], offset[1], offset[2]])
        elif i not in skin and j in skin and (atoms.info['species'][i], atoms.info['species'][j]) in cutoff_min:
            if bondsetting[(atoms.info['species'][i], atoms.info['species'][j])][3]:
                bondpairs_min[i].append([j, offset[0], offset[1], offset[2]])
        elif i in skin and j not in skin and (atoms.info['species'][j], atoms.info['species'][i]) in cutoff_min:
            if bondsetting[(atoms.info['species'][j], atoms.info['species'][i])][3]:
                bondpairs_min[i].append([j, offset[0], offset[1], offset[2]])
    for i, j, offset in zip(nli_max, nlj_max, nlS_max):
        if i not in skin and j not in skin:
            bondpairs_max[i].append([j, offset[0], offset[1], offset[2]])
        elif i not in skin and j in skin and (atoms.info['species'][i], atoms.info['species'][j]) in cutoff_min:
            if bondsetting[(atoms.info['species'][i], atoms.info['species'][j])][3]:
                bondpairs_max[i].append([j, offset[0], offset[1], offset[2]])
        elif i in skin and j not in skin and (atoms.info['species'][j], atoms.info['species'][i]) in cutoff_min:
            if bondsetting[(atoms.info['species'][j], atoms.info['species'][i])][3]:
                bondpairs_max[i].append([j, offset[0], offset[1], offset[2]])
    for i, data in bondpairs_min.items():
        bondpairs = [j for j in bondpairs_max[i] if j not in data]
        bondpairs_max[i] = bondpairs
    bondpairs_max = {i: data  for i, data in bondpairs_max.items() if len(data) > 0}
    ind_atom_skin = list(set(bondpairs_max.keys()) & set(skin))
    print('get_bondpairs 2: {0:10.2f} s'.format(time.time() - tstart))
    if len(ind_atom_skin) == 0:
        return bondpairs_max, atoms_skin
    atoms_skin = atoms[ind_atom_skin]
    atoms_skin.info['species'] = [atoms.info['species'][i] for i in ind_atom_skin]
    return bondpairs_max, atoms_skin

def default_element_prop(element, color_style = "JMOL"):
    """
    """
    from ase.data import chemical_symbols
    from ase.data.colors import jmol_colors, cpk_colors
    from blase.default_data import vesta_color
    element_prop = {}
    number = chemical_symbols.index(element)
    if color_style.upper() == 'JMOL':
        color = jmol_colors[number]
    elif color_style.upper() == 'CPK':
        color = cpk_colors[number]
    elif color_style.upper() == 'VESTA':
        color = vesta_color[element]
    radius = covalent_radii[number]
    element_prop['element'] = element
    element_prop['color'] = color
    element_prop['radius'] = radius
    element_prop['transmit'] = 1.0
    return element_prop

def get_atom_kind(element, positions = [], color_style = "JMOL", scale = [1, 1, 1], props = {}):
    """
    """
    atom_kind = default_element_prop(element, color_style = color_style)
    atom_kind['scale'] = scale
    atom_kind['positions'] = positions
    atom_kind['balltype'] = None
    atom_kind.update(props)
    return atom_kind
def get_bond_kind(element, color_style = "JMOL", props = {}):
    """
    """
    bond_kind = default_element_prop(element, color_style = color_style)
    bond_kind['lengths'] = []
    bond_kind['centers'] = []
    bond_kind['normals'] = []
    bond_kind['verts'] = []
    bond_kind['faces'] = []
    bond_kind.update(props)
    return bond_kind

def get_polyhedra_kind(element, color_style = "JMOL", transmit = 0.8, props = {}):
    polyhedra_kind = default_element_prop(element, color_style = color_style)
    polyhedra_kind['transmit'] = transmit
    vertices = []
    edges = []
    faces = []
    polyhedra_kind.update({'vertices': vertices, 'edges': edges, 'faces': faces})
    lengths = []
    centers = []
    normals = []
    polyhedra_kind['edge_cylinder'] = {'lengths': lengths, 'centers': centers, 'normals': normals}
    polyhedra_kind['edge_cylinder']['color'] = (1.0, 1.0, 1.0)
    polyhedra_kind['edge_cylinder']['transmit'] = 1.0
    polyhedra_kind.update(props)
    return polyhedra_kind

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

def search_boundary(positions, cell, boundary = [[0, 1], [0, 1], [0, 1]], skin = 3.0):
    """
    boundarys: float or list
    """
    from ase.cell import Cell
    from math import floor, ceil
    from time import time
    tstart = time()
    cell = Cell(cell)
    par = cell.cellpar()
    skin = np.array([skin/par[0], skin/par[1], skin/par[2]])
    if isinstance(boundary, float):
        boundary = [[-boundary, 1 + boundary], [-boundary, 1+boundary], [-boundary, 1+boundary]]
    boundary = np.array(boundary)
    boundary[:, 0] += 1e-6
    boundary[:, 1] -= 1e-6
    boundary_skin = boundary.copy()
    boundary_skin[:, 0] -= skin
    boundary_skin[:, 1] += skin
    f = np.floor(boundary_skin)
    c = np.ceil(boundary_skin)
    ib = np.array([f[:, 0], c[:, 1]]).astype(int)
    # print('ib: ', ib)
    M = np.product(ib[1] - ib[0] + 1)
    positions = cell.scaled_positions(positions)
    n = len(positions)
    # print(type(M))
    print('search boundary 1: {0:10.2f} s'.format(time() - tstart))
    tstart = time()
    npositions = np.tile(positions, (M -1,) + (1,) * (len(positions.shape) - 1))
    i0 = 0
    print('search boundary 2: {0:10.2f} s'.format(time() - tstart))
    tstart = time()
    for m0 in range(ib[0, 0], ib[1, 0] + 1):
        for m1 in range(ib[0, 1], ib[1, 1] + 1):
            for m2 in range(ib[0, 2], ib[1, 2] + 1):
                if m0 == 0 and m1 == 0 and m2 == 0: continue
                i1 = i0 + n
                npositions[i0:i1] += (m0, m1, m2)
                i0 = i1
    print('search boundary 3: {0:10.2f} s'.format(time() - tstart))
    tstart = time()
    ind1 = np.array([], dtype=int)
    for i in range(3):
        ind1 = np.append(ind1, np.where(npositions[:, i] < boundary[i][0]))
        ind1 = np.append(ind1, np.where(npositions[:, i] > boundary[i][1]))
    ind2 = np.array([], dtype=int)
    for i in range(3):
        ind2 = np.append(ind2, np.where(npositions[:, i] < boundary_skin[i][0]))
        ind2 = np.append(ind2, np.where(npositions[:, i] > boundary_skin[i][1]))
    print(len(ind1), len(ind2))
    ind1 = list(set(ind1))
    ind2 = list(set(ind2))
    ind3 = [x for x in ind1 if x not in ind2]
    npositions1 = np.delete(npositions, ind1, axis = 0)
    npositions1 = np.dot(npositions1, cell)
    npositions2 = npositions[ind3]
    npositions2 = np.dot(npositions2, cell)
    print('search boundary 4: {0:10.2f} s'.format(time() - tstart))
    return npositions1, npositions2

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
def find_cage(cell, positions, radius, step = 1.0):
    from ase.cell import Cell
    from scipy.spatial import distance

    cell = Cell(cell)
    a, b, c, alpha, beta, gamma = cell.cellpar()
    na, nb, nc = int(a/step),int(b/step), int(c/step)
    x = np.linspace(0, 1, na)
    y = np.linspace(0, 1, nb)
    z = np.linspace(0, 1, nc)
    positions_v = np.vstack(np.meshgrid(x, y, z)).reshape(3,-1).T
    positions_v = np.dot(positions_v, cell)
    dists = distance.cdist(positions_v, positions)
    # dists<3.0
    flag = np.min(dists, axis = 1) >radius
    return positions_v[flag]



if __name__ == "__main__":
    from ase.build import bulk
    from ase.atoms import Atoms
    from ase.io import read
    from ase.visualize import view
    from ase.neighborlist import neighbor_list
    # atoms = bulk('Pt') #, cubic = True)
    # atoms = bulk('Pt', cubic = True)
    atoms = read('docs/source/_static/datas/tio2.cif')
    nli_min, nlj_min, nlS_min = neighbor_list('ijS', atoms, cutoff={('Ti', 'O'): 2.5, ('O', 'O'): 1.4}, self_interaction=False)
    print(nli_min)
    print(nlj_min)
    view(atoms)
    # positions = find_cage(atoms.cell, atoms.positions, 9.0, step = 1.0)
    # vatoms = Atoms(['Au']*len(positions), positions=positions)
    positions, positions2 = search_boundary(atoms.positions, atoms.cell, boundary=[[-0.01, 1.01], [-0.01, 1.01], [-0.01, 1.01]])
    vatoms = Atoms(['Au']*len(positions), positions=positions)
    vatoms2 = Atoms(['Au']*len(positions2), positions=positions2)
    # print(len(vatoms2))
    atoms = atoms + vatoms + vatoms2
    atoms.pbc = False
    # print(vatoms)
    # view([atoms, vatoms, vatoms2])
    # view(vatoms2)
    # atoms.pbc = False

