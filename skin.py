import numpy as np
from math import floor, ceil
from time import time
from ase import Atom, Atoms
from ase.cell import Cell



def build_core(ib0, n):
    #--------
    # core
    npositions = []
    ind = list(range(n))
    for m0 in range(ib0[0, 0], ib0[1, 0]):
        for m1 in range(ib0[0, 1], ib0[1, 1]):
            for m2 in range(ib0[0, 2], ib0[1, 2]):
                if m0 == 0 and m1 == 0 and m2 == 0: continue
                npositions.append([ind, [m0, m1, m2]])
    print('core: ', len(npositions))
    return npositions

def build_facet(ib0, ib1, ind_bd_cells):
    # -x, facet
    npositions = []
    dcell = np.array([0, 0, 0])
    for x, yz in {0: [1, 2], 1:[0, 2], 2:[0, 1]}.items():
        for i in range(2):
            ind = ind_bd_cells[i][x]
            if len(ind) == 0: continue
            if i == 0:
                m0 = ib1[i, x]
            else:
                m0 = ib1[i, x] - 1
            for m1 in range(ib0[0, yz[0]], ib0[1, yz[0]]):
                for m2 in range(ib0[0, yz[1]], ib0[1, yz[1]]):
                    dcell[x] = m0
                    dcell[yz[0]] = m1
                    dcell[yz[1]] = m2
                    print('dcell: ', dcell)
                    npositions.append([ind, dcell.copy()])
                    print('x: ', x, 'i: ', i, 'ind: ', ind)
    print('facet: ', len(npositions))
    return npositions
def build_edge(ib0, ib1, ind_bd_cells):
    """
    """
    # -xy, edge
    npositions = []
    dcell = np.array([0, 0, 0], dtype=int)
    for xy, z in {(0, 1): 2, (0, 2):1, (1, 2):0}.items():
        for i in range(2):
            for j in range(2):
                ind = list(set(ind_bd_cells[i][xy[0]]) & set(ind_bd_cells[j][xy[1]]))
                if len(ind) == 0: continue
                m0 = ib1[i, xy[0]] - i
                m1 = ib1[j, xy[1]] - j
                # print('ind: ', ind)
                for m2 in range(ib0[0, z], ib0[1, z]):
                    dcell[xy[0]] = m0
                    dcell[xy[1]] = m1
                    dcell[z] = m2
                    print('dcell: ', dcell)
                    npositions.append([ind, dcell.copy()])
                    print('xy: ', xy, 'i, j:', i, j, 'ind: ', ind)
    print('edge: ', len(npositions))
    return npositions
def build_corner(ib0, ib1, ind_bd_cells):
    """
    """
    # -xyz, corner
    m = ib0[1] - ib0[0]
    npositions = []
    dcell = np.array([0, 0, 0], dtype=int)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                ind = list(set(ind_bd_cells[i][0]) & set(ind_bd_cells[j][1]) & set(ind_bd_cells[k][2]))
                if len(ind) == 0: continue
                m0 = ib1[i, 0] - i
                m1 = ib1[j, 1] - j
                m2 = ib1[k, 2] - k
                # print('ind: ', ind)
                dcell[0] = m0
                dcell[1] = m1
                dcell[2] = m2
                # print('dcell: ', dcell)
                npositions.append([ind, dcell.copy()])
                # print('i, j:', i, j, 'ind: ', ind, 'M1: ', M1, npositions)
    print('corner: ', len(npositions))
    return npositions

def search_boundary_cell(positions, boundary):
    """
    find boundary atoms inside unit cell
    """
    f0 = np.ceil(boundary)
    c0 = np.floor(boundary)
    f1 = np.floor(boundary)
    c1 = np.ceil(boundary)
    ib0 = np.array([f0[:, 0], c0[:, 1]]).astype(int)
    ib1 = np.array([f1[:, 0], c1[:, 1]]).astype(int)
    print('ib0: ', ib0)
    print('ib1: ', ib1)
    ind_bd_cells = [[0, 0, 0], [0, 0, 0]]
    c = 1 - (np.ceil(boundary) - boundary)  # 1 - (0 - (-0.6)) = 0.4,   > 0.4, < 1
    f = boundary - np.floor(boundary)  # 1.6 - 1 = 0.6, > 0, < 0.6, a > 0
    dboundary = [0, 0]
    for i in range(3):
        # boundary < 0 or >0
        for j in range(2):
            if boundary[i][j] <= 0:
                dboundary[j] = 1 - (np.ceil(boundary[i][j]) - boundary[i][j])
            else:
                dboundary[j] = boundary[i][j] - np.floor(boundary[i][j]) 
        print('dboundary: ', dboundary)
        # boundary in same side or not
        if ib0[0][i] - ib0[1][i] == 1:
            # same side, both condition
            if boundary[i][j] > 0:
                ind = list(np.where(np.logical_and(positions[:, i] > dboundary[0], positions[:, i] < dboundary[1]))[0])
                ind_bd_cells[0][i] = []
                ind_bd_cells[1][i] = ind
            else:
                ind = list(np.where(np.logical_and(positions[:, i] > dboundary[0], positions[:, i] < dboundary[1]))[0])
                ind_bd_cells[0][i] = ind
                ind_bd_cells[1][i] = []
        else:
            ind = list(np.where(np.logical_and(positions[:, i] > dboundary[0], positions[:, i] <= 1))[0])
            ind_bd_cells[0][i] = ind
            ind = list(np.where(np.logical_and(positions[:, i] < dboundary[1], positions[:, i] >= 0))[0])
            ind_bd_cells[1][i] = ind
            
    print('ind_bd_cells: ', ind_bd_cells)
    # both condition should be fullfilled
    for xyz in [[0, 1, 2], [1, 0, 2], [2, 0, 1]]:
        x = xyz[0]
        for yz in xyz[1:]:
            if ib0[0][yz] - ib0[1][yz] == 1:
                if boundary[yz][0] > 0:
                    print('x, yz: ', x, yz, '  > 0')
                    ind_bd_cells[0][x] = set(ind_bd_cells[0][x]) & set(ind_bd_cells[1][yz])
                    ind_bd_cells[1][x] = set(ind_bd_cells[1][x]) & set(ind_bd_cells[1][yz])
                else:
                    print('x, yz: ', x, yz, '  < 0')
                    ind_bd_cells[0][x] = set(ind_bd_cells[0][x]) & set(ind_bd_cells[0][yz])
                    ind_bd_cells[1][x] = set(ind_bd_cells[1][x]) & set(ind_bd_cells[0][yz])

        ind_bd_cells[0][x] = list(ind_bd_cells[0][x] )
        ind_bd_cells[1][x] = list(ind_bd_cells[1][x] )
    print('ind_bd_cells: ', ind_bd_cells)
    return ib0, ib1, ind_bd_cells
    
def search_boundary(positions, cell, boundary):
    """
    boundarys: float or list

    """
    from time import time
    cell = Cell.new(cell)
    tstart0 = time()
    if isinstance(boundary, float):
        boundary = np.array([[-boundary, 0], [1, 1 + boundary]], [[-boundary, 0], [1, 1+boundary]], [[-boundary, 0], [1, 1+boundary]])
    boundary = np.array(boundary)
    print('boundary: ', boundary)
    positions = cell.scaled_positions(positions)
    ib0, ib1, ind_bd_cells = search_boundary_cell(positions, boundary)
    npositions0 = build_core(ib0, len(positions))
    #--------
    #--------
    npositions1 = build_facet(ib0, ib1, ind_bd_cells)
    #----------------------------------------------------
    npositions2 = build_edge(ib0, ib1, ind_bd_cells)
    #----------------------------------------------------
    npositions3 = build_corner(ib0, ib1, ind_bd_cells)
    #------------------------------------------------------------
    npositions0.extend(npositions1)
    npositions0.extend(npositions2)
    npositions0.extend(npositions3)
    print('search boundary: {0:10.2f} s'.format(time() - tstart0))
    return npositions0
def index2positions(index, positions, cell):
    """
    """
    cell = Cell.new(cell)
    positions = cell.scaled_positions(positions)
    npositions = []
    for data in index:
        ind, offset = data
        temp = positions[ind]
        temp += offset
        npositions.extend(temp)
    print('npositions: ', len(npositions))
    npositions = np.dot(npositions, cell)
    return npositions
def search_skin(positions, cell, boundary, skin = 2.0):
    """
    skins: float or list

    """
    from time import time
    tstart0 = time()
    boundary0 = np.array(boundary)
    print('boundary 0: ', boundary0)
    cell = Cell.new(cell)
    par = cell.cellpar()
    skin = np.array([skin/par[i] for i in range(3)])
    print('skin: ', skin)
    positions = cell.scaled_positions(positions)
    #--------
    # facet
    npositions = []
    for i in range(3):
        print('skin: ', i)
        for j in range(2):
            print('  skin: ', j)
            boundary = boundary0.copy()
            if j == 0:
                boundary[i][0] = boundary0[i, 0] - skin[i]
                boundary[i][1] = boundary0[i, 0]
            if j == 1:
                boundary[i][0] = boundary0[i, 1]
                boundary[i][1] = boundary0[i, 1] + skin[i]
            temp = search_boundary(positions, cell, boundary)
            npositions.extend(temp)
    #--------
    # edge, six
    for i in range(3):
        for j in range(3):
            for k in range(2):
                boundary = boundary0.copy()
                if j == 0:
                    boundary[i][0] = boundary0[i, 0] - skin[i]
                    boundary[i][1] = boundary0[i, 0]
                if j == 1:
                    boundary[i][0] = boundary0[i, 1]
                    boundary[i][1] = boundary0[i, 1] + skin[i]
                temp = search_boundary(positions, cell, boundary)
                # npositions.extend(temp)
    #----------------------------------------------------
    # corner
    #----------------------------------------------------
    #------------------------------------------------------------
    print('search skin: {0:10.2f} s'.format(time() - tstart0))
    return npositions

class Boundary:
    """
    Boundary object.

    Example:
      cl = Boundary(atoms, d, index)
      cl.build()
    """

    def __init__(self, atoms, boundary_list, rotate_atoms = False):
        self.atoms = atoms
        self.natoms = len(atoms)
        self.cell = atoms.cell
        self.boundary_list = boundary_list
        # self.d = d
        # self.index = index
        self.rotate_atoms = rotate_atoms
        #
    def build(self, ):
        for boundary in self.boundary_list:
            print(boundary)
            self.cut(**boundary)
        return self.atoms
    def cut(self, atoms = None, d = None, index = None, direction = 1):
        """
        """
        if not atoms:
            atoms = self.atoms
        cell = atoms.cell
        normal = self.get_plane(d, index, cell)
        print(normal, d)
        # a*x + b*y + c*z - d = 0
        mask = []
        natoms = len(atoms)
        for i in range(natoms):
            v = atoms[i].position.dot(normal) - d
            # print(v)
            if v*direction > 0:
                mask.append(i)
        # view(atoms[mask])
        del atoms[mask]
        self.atoms = atoms
        if self.rotate_atoms:
            atoms = self.rotate()
    def rotate(self, atoms = None, index = None):
        """
        rotate normal of plane to z axis
        """
        import scipy
        if not atoms:
            atoms = self.atoms
        if not index:
            index = self.index
        cell = atoms.cell
        normal, d = self.get_plane(self.d, self.index, cell)
        vec = np.cross([0.0000014159, 0.000001951, 1], index)
        vec = vec/np.linalg.norm(vec)
        ang = np.arccos(normal[2])
        vec = ang*vec
        r = scipy.spatial.transform.Rotation.from_rotvec(vec)
        if scipy.version.version >= '1.4':
            mat = r.as_matrix()
        else: 
            mat = r.as_dcm()
        # print(mat)
        atoms.positions = atoms.positions.dot(mat)
        # atoms.cell = atoms.cell.dot(mat)
        self.atoms = atoms
        return atoms
    def get_plane(self, d, index = None, cell = None):
        '''
        plane equation: three point and distance from origin
        return normal and distance
        # a*x + b*y + c*z - d = 0
        '''
        index = [1.0/(index[i] + 0.000001) for i in range(3)]
        index = np.array(index)
        index = index/np.linalg.norm(index)
        points = cell*index
        # print(points)
        # x 
        v1 = points[1] - points[0]
        v2 = points[2] - points[0]
        normal = np.cross(v1, v2)
        print(v1, v2, normal)
        normal = normal/np.linalg.norm(normal)
        a, b, c = normal
        # d = np.dot(normal, points[2])
        return normal

if __name__ == "__main__":
    from ase.io import read, write
    from ase.build import molecule, bulk
    from ase.visualize import view
    # atoms = bulk('Pt', cubic = True)
    # atoms = atoms*[6, 6, 6]
    # cl = Boundary(atoms, boundary_list = [{'d': 10.0, 'index': [2, 2, 1]}])
    # atoms = cl.build()
    # view(atoms)
    atoms = bulk('Pt', cubic = True)
    atoms.write('pt.in')
    # atoms = atoms*[2, 2, 1]
    # atoms = read('docs/source/_static/datas/tio2.cif')
    # atoms = read('docs/source/_static/datas/mof-5.cif')
    # atoms.positions -= atoms.get_center_of_mass()
    # index = search_boundary(atoms.positions, atoms.cell, boundary=[[-0.6, 1.6], [-0.6, 1.6], [0, 1]])
    index = search_skin(atoms.positions, atoms.cell, boundary=[[-0.6, 1.6], [-0.6, 1.6], [0, 1]])
    positions = index2positions(index, atoms.positions, atoms.cell)
    vatoms = Atoms(['Au']*len(positions), positions=positions)
    # vatoms2 = Atoms(['Au']*len(positions2), positions=positions2)
    # atoms = atoms*[3, 3, 3]
    # atoms.positions -= atoms.get_center_of_mass()
    atoms = atoms + vatoms
    atoms.pbc = False
    print(vatoms)
    view(atoms)
    # view([atoms, vatoms, vatoms2])
    # atoms.pbc = False

