import numpy as np
from ase import Atoms, Atom, atom
from ase.data import covalent_radii, atomic_numbers, chemical_symbols
from ase.visualize import view
import time
from ase.neighborlist import neighbor_list

    
def default_element_prop(element, color_style = "JMOL"):
    """
    """
    from ase.data import chemical_symbols
    from ase.data.colors import jmol_colors, cpk_colors
    from default_data import vesta_color
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
    atoms = bulk('Pt', cubic = True)
    # atoms = atoms*[2, 2, 1]
    # atoms = read('docs/source/_static/datas/tio2.cif')
    # atoms = read('docs/source/_static/datas/mof-5.cif')
    # nli_min, nlj_min, nlS_min = neighbor_list('ijS', atoms, cutoff={('Ti', 'O'): 2.5, ('O', 'O'): 1.4}, self_interaction=False)
    # print(nli_min)
    # print(nlj_min)
    # view(atoms)
    # positions = find_cage(atoms.cell, atoms.positions, 9.0, step = 1.0)
    # vatoms = Atoms(['Au']*len(positions), positions=positions)
    