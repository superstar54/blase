import numpy as np
from time import time

def search_boundary(positions, cell, boundary = [[0, 1], [0, 1], [0, 1]], skin = 3.0):
    """
    boundarys: float or list
    """
    from ase.cell import Cell
    from math import floor, ceil
    tstart = time()
    cell = Cell(cell)
    par = cell.cellpar()
    skin = np.array([skin/par[0], skin/par[1], skin/par[2]])
    if isinstance(boundary, float):
        boundary = [[-boundary, 1 + boundary], [-boundary, 1+boundary], [-boundary, 1+boundary]]
    boundary = np.array(boundary)
    boundary_skin = boundary.copy()
    boundary_skin[:, 0] -= skin
    boundary_skin[:, 1] += skin
    #
    boundary_i = boundary.copy()
    boundary_i[:, 0] -= skin
    boundary_i[:, 1] += skin
    f = np.floor(boundary)
    c = np.ceil(boundary)
    ib = np.array([f[:, 0], c[:, 1]]).astype(int)
    M = np.product(ib[1] - ib[0] + 1)
    positions = cell.scaled_positions(positions)
    n = len(positions)
    npositions = np.tile(positions, (M - 1,) + (1,) * (len(positions.shape) - 1))
    i0 = 0
    # repeat the positions so that
    # it completely covers the boundary
    for m0 in range(ib[0, 0], ib[1, 0] + 1):
        for m1 in range(ib[0, 1], ib[1, 1] + 1):
            for m2 in range(ib[0, 2], ib[1, 2] + 1):
                if m0 == 0 and m1 == 0 and m2 == 0: continue
                i1 = i0 + n
                npositions[i0:i1] += (m0, m1, m2)
                i0 = i1
    # boundary
    ind1 =  np.where((npositions[:, 0] > boundary[0][0]) & (npositions[:, 0] < boundary[0][1]) \
                & (npositions[:, 1] > boundary[1][0]) & (npositions[:, 1] < boundary[1][1]) \
                & (npositions[:, 2] > boundary[2][0]) & (npositions[:, 2] < boundary[2][1]))
    #
    # boundary
    ind2 =  np.where((npositions[:, 0] > boundary_skin[0][0]) & (npositions[:, 0] < boundary_skin[0][1]) \
                & (npositions[:, 1] > boundary_skin[1][0]) & (npositions[:, 1] < boundary_skin[1][1])  #\
                & (npositions[:, 2] > boundary_skin[2][0]) & (npositions[:, 2] < boundary_skin[2][1]))
    # boundary
    ind3 =  np.where((npositions[:, 0] > boundary_skin[0][0]) & (npositions[:, 0] < boundary_skin[0][1]) \
                & (npositions[:, 1] > boundary_skin[1][0]) & (npositions[:, 1] < boundary_skin[1][1])  #\
                & (npositions[:, 2] > boundary_skin[2][0]) & (npositions[:, 2] < boundary_skin[2][1]))
    ind1 = list(set(ind1[0]))
    ind2 = list(set(ind2[0]))
    ind2 = list(set(ind2) - set(ind1))
    npositions1 = npositions[ind1]
    npositions1 = np.dot(npositions1, cell)
    npositions2 = npositions[ind2]
    npositions2 = np.dot(npositions2, cell)
    # print('search boundary: {0:10.2f} s'.format(time() - tstart))
    return npositions1, npositions2

def search_skin(atoms, bondsetting):
    from blase.bondsetting import build_bondlists
    from ase import Atoms
    cutoff = {}
    for b in bondsetting:
        if b.search > 0:
            cutoff[(b.symbol1, '%s_skin'%b.symbol2)] = [b.min, b.max]
    bondlist = build_bondlists(atoms, cutoff)
    if len(bondlist) == 0: return Atoms()
    index = bondlist[:, 1]
    atoms_skin = atoms[index]
    atoms_skin.info['species'] = np.array(atoms.info['species'])[index]
    return atoms_skin


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
    from ase import Atom, Atoms
    # atoms = bulk('Pt', cubic = True)
    # atoms = atoms*[6, 6, 6]
    # cl = Boundary(atoms, boundary_list = [{'d': 10.0, 'index': [2, 2, 1]}])
    # atoms = cl.build()
    # view(atoms)
    atoms = bulk('Pt', cubic = True)
    atoms.write('pt.in')
    # view(atoms)
    # atoms = atoms*[3, 3, 3]
    # atoms = read('docs/source/_static/datas/tio2.cif')
    # atoms = read('docs/source/_static/datas/mof-5.cif')
    # atoms.positions -= atoms.get_center_of_mass()
    positions1, positions2 = search_boundary(atoms.positions, atoms.cell, boundary=[[-0.6, 1.6], [-0.6, 1.6], [-0.6, 1.6]])
    skin = Atoms('Au'*len(positions2), positions = positions2)
    view(atoms + skin)
    # bondsetting = {('Ti', 'O'): [0, 2.5, False, False], ('O', 'O'): [0, 1.5, False, False]}
    # bondpairs_boundary = search_skin(atoms, bondsetting, boundary=[[-0.6, 1.6], [-0.6, 1.6], [-0.6, 1.6]])
    # print(bondpairs_boundary)
