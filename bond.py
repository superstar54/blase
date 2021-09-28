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
    atoms = molecule('NH3')
    atoms[1].z -= 0.5
    bondsetting = {('H', 'N'): [0, 1.5, False, False]}
    bondlists = get_bondlists(atoms, bondsetting)
    datas = calc_bond_data(atoms, bondlists)
    source = ([[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]], [[0, 1, 2, 3]])
    cylinder_mesh_from_instance_vec(datas['N']['centers'], datas['N']['normals'], datas['N']['lengths'], 0.1, source)
    # print(bondlists)
    