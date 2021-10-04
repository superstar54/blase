import numpy as np
from time import time
from pprint import pprint


def get_bondtable(specieslist, cutoff = 1.2, add_bonds = {}, remove_bonds = {}, polyhedra_dict = {}):
    """
    """
    from blase.default_data import default_bonds
    from ase.data import chemical_symbols, covalent_radii
    bondtable = {}
    for species1 in specieslist:
        search = False
        polyhedra = False
        radius1 = cutoff * covalent_radii[chemical_symbols.index(species1.split('_')[0])]
        for species2 in specieslist:
            if species2 not in default_bonds[species1.split('_')[0]]: continue
            radius2 = cutoff * covalent_radii[chemical_symbols.index(species2.split('_')[0])]
            bondmax = radius1 + radius2
            if species1 in polyhedra_dict:
                if species2 in polyhedra_dict[species1]:
                    polyhedra = True
            bondtable[(species1, species2)] = [0.5, bondmax, polyhedra, search]
    bondtable.update(add_bonds)
    for key in remove_bonds:
        bondtable.pop(key)
    return bondtable

def build_bondlists(atoms, bondsetting):
    """
    The default bonds are stored in 'default_bonds'
    Get all pairs of bonding atoms
    remove_bonds
    """
    from blase.neighborlist import neighbor_list
    if len(bondsetting) == 0: return {}
    cutoff_min = {}
    cutoff = {}
    for key, data in bondsetting.items():
        cutoff_min[key] = data[0]
        cutoff[key] = data[1]
    #
    tstart = time()
    nli, nlj, nlS = neighbor_list('ijS', atoms, cutoff, cutoff_min, self_interaction=False)
    bondlists = np.append(np.array([nli, nlj], dtype=int).T, np.array(nlS, dtype=int), axis = 1)
    # print('build_bondlists: {0:10.2f} s'.format(time() - tstart))
    return bondlists

def build_polyhedralists(atoms, bondlists, bondsetting, color_style = "JMOL", transmit = 0.8):
        """
        Two modes:
        (1) Search atoms bonded to kind
        polyhedra_dict: {'kind': ligands}
        """
        from scipy.spatial import ConvexHull
        from blase.tools import get_polyhedra_kind
        tstart = time()
        if 'species' not in atoms.info:
            atoms.info['species'] = atoms.get_chemical_symbols()
        speciesarray = np.array(atoms.info['species'])
        speciesarray0 = speciesarray[bondlists[:, 0]]
        speciesarray1 = speciesarray[bondlists[:, 1]]
        positions = atoms.positions
        polyhedra_kinds = {}
        polyhedra_dict = {}
        ind0 = []
        for (spi, spj), data in bondsetting.items():
            if data[2]:
                if spi not in polyhedra_dict: polyhedra_dict[spi] = []
                polyhedra_dict[spi].append(spj)
                ind0.append(list(np.where((speciesarray0 == spi) & (speciesarray1 == spj))[0]))
        # loop center atoms
        npositions = positions[bondlists[:, 1]] + np.dot(bondlists[:, 2:5], atoms.cell)
        # data[2] is ture
        for kind, ligand in polyhedra_dict.items():
            # print(kind, ligand)
            if kind not in polyhedra_kinds.keys():
                element = kind.split('_')[0]
                polyhedra_kinds[kind] = get_polyhedra_kind(element, color_style=color_style)
            inds = np.where(speciesarray == kind)[0]
            for ind in inds:
                vertice = npositions[bondlists[:, 0] == ind]
                nverts = len(vertice)
                if nverts >3:
                    # search convex polyhedra
                    hull = ConvexHull(vertice)
                    face = hull.simplices
                    #
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
            # print(polyhedra_kinds)
        # print('build_polyhedralists: {0:10.2f} s'.format(time() - tstart))
        return polyhedra_kinds

def search_skin(atoms, bondsetting, bondlists, skin = []):
    """
    The default bonds are stored in 'default_bonds'
    Get all pairs of bonding atoms
    remove_bonds
    """
    from ase import Atoms
    tstart = time()
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
    # print('search skin : {0:10.2f} s'.format(time() - tstart))
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
    atoms = atoms*[4, 4, 4]
    # positions1, offsets1, positions2, offsets2 = search_boundary(atoms.positions, atoms.cell, boundary=[[-0.6, 1.6], [-0.6, 1.6], [-0.6, 1.6]])
    bondsetting = {('Ti', 'O'): [0, 2.5, True, False], ('O', 'O'): [0, 1.5, False, False]}
#     bondsetting = {('H' ,   'O' ):    [0.500,   1.164,      False,        False      ],
# ('O' ,   'H' ):    [0.500,   1.164,      False,        False      ],
# ('O' ,   'O' ):    [0.500,   1.584,      False,        False      ],
# ('Zn',    'H'):     [0.500,   1.836,      False,        False      ],
# ('Zn',    'O'):     [0.000,   2.500,      True ,        False      ],
# ('C' ,   'H' ):    [0.000,   1.400,      False,        False      ],
# ('C' ,   'O' ):    [0.500,   1.704,      False,        False      ],
# ('C' ,   'C' ):    [0.500,   1.824,      False,        False],}
    bondlists = build_bondlists(atoms, bondsetting)
    datas = build_polyhedralists(atoms, bondlists, bondsetting)
    