import numpy as np
from ase import Atoms, Atom
from ase.data import covalent_radii, atomic_numbers, chemical_symbols
from ase.visualize import view
import time

def get_bondpairs(atoms, bondsetting):
    """
    The default bonds are stored in 'default_bonds'
    Get all pairs of bonding atoms
    remove_bonds
    """
    from ase.neighborlist import neighbor_list
    tstart = time.time()
    cutoff = {}
    for key, data in bondsetting.items():
        cutoff[key] = data[0]
    nli, nlj, nlS = neighbor_list('ijS', atoms, cutoff=cutoff, self_interaction=False)
    bondpairs = {i: [] for i in set(nli)}
    if 'species' not in atoms.info:
        atoms.info['species'] = atoms.get_chemical_symbols()
    for i, j, offset in zip(nli, nlj, nlS):
        bondpairs[i].append([j, offset])
    print('get_bondpairs: {0:10.2f} s'.format(time.time() - tstart))
    return bondpairs

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

def get_atom_kind(element, positions, color_style = "JMOL", scale = [1, 1, 1], props = {}):
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
    pass