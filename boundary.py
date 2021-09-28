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
    # print('boundary: ', boundary)
    # print('boundary_skin: ', boundary_skin)
    f = np.floor(boundary)
    c = np.ceil(boundary)
    ib = np.array([f[:, 0], c[:, 1]]).astype(int)
    M = np.product(ib[1] - ib[0] + 1)
    positions = cell.scaled_positions(positions)
    n = len(positions)
    npositions = np.tile(positions, (M - 1,) + (1,) * (len(positions.shape) - 1))
    i0 = 0
    offsets = np.zeros((M*n, 4), dtype=int)
    ind0 = np.arange(n).reshape(-1, 1)
    #
    # print('ib: ', ib)
    for m0 in range(ib[0, 0], ib[1, 0] + 1):
        for m1 in range(ib[0, 1], ib[1, 1] + 1):
            for m2 in range(ib[0, 2], ib[1, 2] + 1):
                if m0 == 0 and m1 == 0 and m2 == 0: continue
                i1 = i0 + n
                npositions[i0:i1] += (m0, m1, m2)
                offsets[i0:i1] = np.append(ind0, np.array([[m0, m1, m2]]*n), axis = 1)
                i0 = i1
    #
    ind1 =  np.where((npositions[:, 0] > boundary[0][0]) & (npositions[:, 0] < boundary[0][1]) \
                & (npositions[:, 1] > boundary[1][0]) & (npositions[:, 1] < boundary[1][1]) \
                & (npositions[:, 2] > boundary[2][0]) & (npositions[:, 2] < boundary[2][1]))
    #
    # print(npositions)
    npositions = np.append(npositions, positions, axis = 0)
    offsets[i0:i0+n] = np.append(ind0, np.array([[m0, m1, m2]]*n), axis = 1)
    ind2 =  np.where((npositions[:, 0] > boundary_skin[0][0]) & (npositions[:, 0] < boundary_skin[0][1]) \
                & (npositions[:, 1] > boundary_skin[1][0]) & (npositions[:, 1] < boundary_skin[1][1])  #\
                & (npositions[:, 2] > boundary_skin[2][0]) & (npositions[:, 2] < boundary_skin[2][1]))
    ind1 = list(set(ind1[0]))
    ind2 = list(set(ind2[0]))
    ind2 = list(set(ind2) - set(ind1))
    # print(ind1, ind2)
    npositions1 = npositions[ind1]
    offsets1 = offsets[ind1]
    npositions1 = np.dot(npositions1, cell)
    npositions2 = npositions[ind2]
    offsets2 = offsets[ind2]
    npositions2 = np.dot(npositions2, cell)
    print('search boundary: {0:10.2f} s'.format(time() - tstart))
    return npositions1, offsets1, npositions2, offsets2

def search_skin(atoms, bondsetting, boundary):
    """
    build bond
    Then search bond
    """
    from bond import get_bondpairs, get_bondpairs_boundary
    # get bondpairs from primaril cell
    bondpairs = get_bondpairs(atoms, bondsetting)
    # print(bondpairs)
    # search boundary
    npositions1, offsets1, npositions2, offsets2 = search_boundary(atoms.positions, atoms.cell, boundary=boundary)
    # print(offsets1)
    # build bond paris for all atoms
    bondpairs_boundary = get_bondpairs_boundary(bondpairs, offsets2)
    # check search bond
    
    return bondpairs_boundary
    



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
    positions1, offsets1, positions2, offsets2 = search_boundary(atoms.positions, atoms.cell, boundary=[[-0.6, 1.6], [-0.6, 1.6], [-0.6, 1.6]])
    skin = Atoms('Au'*len(positions2), positions = positions2)
    view(atoms + skin)
    # bondsetting = {('Ti', 'O'): [0, 2.5, False, False], ('O', 'O'): [0, 1.5, False, False]}
    # bondpairs_boundary = search_skin(atoms, bondsetting, boundary=[[-0.6, 1.6], [-0.6, 1.6], [-0.6, 1.6]])
    # print(bondpairs_boundary)
