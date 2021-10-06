"""
"""
import bpy
from blase.butils import object_mode
import numpy as np
from time import time
from pprint import pprint

class Bondsetting():
    """
    Bondsetting object

    The Bondsetting object store the bondpair information.

    Parameters:

    batoms: Batoms object
        The batoms object that a Bondsetting belong to.
    cutoff: float
        Cutoff used to calculate the maxmium bondlength for bond pairs.
    """
    def __init__(self, label, cutoff = 1.3, color_style = 'JMOL') -> None:
        self.label = label
        self.cutoff = cutoff
        self.color_style = color_style
        self.set_default(self.species, cutoff)
    def get_data(self):
        data = {}
        coll = bpy.data.collections[self.label]
        for b in coll.bond:
            data[(b.symbol1, b.symbol2)] = [b.min, b.max, b.search, b.polyhedra, 
                                            np.array(b.color1), np.array(b.color2), b.bondlinewidth]
        return data
    def get_color(self):
        color = {}
        coll = bpy.data.collections[self.label]
        for b in coll.bond:
            color[(b.symbol1, b.symbol2)] = [np.array(b.color1), np.array(b.color2)]
        return color
    @property
    def species(self):
        return self.get_species()
    def get_species(self):
        """
        read species from collection.
        """
        species = []
        coll_atom = bpy.data.collections['%s_atom'%self.label]
        for ba in coll_atom.objects:
            species.append(ba.species)
        return species
    @property
    def data(self):
        return self.get_data()
    @property
    def color(self):
        return self.get_color()
    @color.setter
    def color(self, color):
        self.set_color(color)
    def set_color(self, color):
        bond, color = color
        data = self.data[bond]
        data[4] = color[0]
        data[5] = color[1]
        self[bond] = data
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Bondpair      min     max   Search_bond    Polyhedra \n'
        data = self.data
        for key, value in data.items():
            s += '{0:5s} {1:5s} {2:4.3f}   {3:4.3f}      {4:10s}   {5:10s} \n'.format(\
                key[0], key[1], value[0], value[1], str(value[2]), str(value[3]),)
        s += '-'*60 + '\n'
        return s
    def __getitem__(self, index):
        return self.data[index]
    def __setitem__(self, index, value):
        """
        Add bondpair one by one
        """
        coll = bpy.data.collections[self.label]
        flag = False
        for b in coll.bond:
            if (b.symbol1, b.symbol2) == index:
                b.min = value[0]
                b.max = value[1]
                b.search = value[2]
                b.polyhedra = value[3]
                if len(value) == 7:
                    b.color1 = value[4]
                    b.color2 = value[5]
                    b.bondlinewidth = value[6]
                flag = True
        if not flag:
            bond = coll.bond.add()
            bond.symbol1 = index[0]
            bond.symbol2 = index[1]
            bond.min = value[0]
            bond.max = value[1]
            bond.search = value[2]
            bond.polyhedra = value[3]
            if len(value) == 7:
                bond.color1 = value[4]
                bond.color2 = value[5]
                bond.bondlinewidth = value[6]
    def copy(self, label):
        object_mode()
        bondsetting = Bondsetting(label)
        for key, value in self.data.items():
            bondsetting[key] = value
        return bondsetting
    def set_default(self, species, cutoff = 1.3):
        """
        """
        bondtable = get_bondtable(species, cutoff=cutoff)
        for key, value in bondtable.items():
            self[key] = value
    def __add__(self, other):
        self += other
        return self
    def __iadd__(self, other):
        self.extend(other)
        return self
    def extend(self, other):
        for key, value in other.data.items():
            self[key] = value
        # new
        species1 = set(self.species)
        species2 = set(other.species)
        nspecies1 = species1 - species2
        nspecies2 = species2 - species1
        for sp1 in nspecies1:
            for sp2 in nspecies2:
                self.set_default([sp1, sp2], self.cutoff)
    def add_bonds(self, bondpair):
        for key in bondpair:
            bondtable.pop(key)
    def remove_bonds(self, bondpair):
        for key in bondpair:
            bondtable.pop(key)
    def hydrogen_bond(self, ):
        if self.hydrogen_bond:
            self.hydrogen_bondlist = build_bondlists(self.atoms, cutoff = {('O', 'H'): self.hydrogen_bond})
        pass
    def add_bondlist(self, ):
        pass
    def add_polyhedra_dict(self, ):
        pass
    
    

def get_bondtable(specieslist, cutoff = 1.3, color_style = "JMOL"):
    """
    """
    from blase.default_data import default_bonds
    from ase.data import chemical_symbols, covalent_radii
    from blase.tools import default_element_prop
    bondtable = {}
    for species1 in specieslist:
        element1 = species1.split('_')[0]
        prop1 = default_element_prop(element1, color_style=color_style)
        color1 = np.append(prop1['color'], [1.0])
        radius1 = cutoff * prop1['radius']
        for species2 in specieslist:
            element2 = species2.split('_')[0]
            pair = (element1, element2)
            if pair not in default_bonds: continue
            prop2 = default_element_prop(element2, color_style=color_style)
            color2 = np.append(prop2['color'], [1.0])
            radius2 = cutoff * prop2['radius']
            bondmax = radius1 + radius2
            bondtable[(species1, species2)] = [0.5, bondmax, default_bonds[pair][0], default_bonds[pair][1], color1, color2, 0.10]
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
    for (spi, spj), data in bondsetting.items():
        eli = spi.split('_')[0]
        elj = spj.split('_')[0]
        key = (eli, elj)
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
        speciesarray_bondlist = speciesarray[bondlists[:, 0]]
        positions = atoms.positions
        polyhedra_kinds = {}
        polyhedra_dict = {}
        for (spi, spj), data in bondsetting.items():
            if data[3]:
                if spi not in polyhedra_dict: polyhedra_dict[spi] = []
                polyhedra_dict[spi].append(spj)
        # loop center atoms
        npositions = positions[bondlists[:, 1]] + np.dot(bondlists[:, 2:5], atoms.cell)
        for spi, spjs in polyhedra_dict.items():
            element = spi.split('_')[0]
            polyhedra_kind = get_polyhedra_kind(element, color_style=color_style)
            indis = np.where(speciesarray == spi)[0]
            for indi in indis:
                vertices = []
                indjs = np.where(bondlists[:, 0] == indi)[0]
                for indj in indjs:
                    if speciesarray[bondlists[indj, 1]] in spjs:
                        vertices.append(npositions[indj])
                nverts = len(vertices)
                if nverts >= 4:
                    # search convex polyhedra
                    hull = ConvexHull(vertices)
                    face = hull.simplices
                    nverts = len(polyhedra_kind['vertices'])
                    face = face + nverts
                    edge = []
                    for f in face:
                        edge.append([f[0], f[1]])
                        edge.append([f[0], f[2]])
                        edge.append([f[1], f[2]])
                    polyhedra_kind['vertices'] = polyhedra_kind['vertices'] + list(vertices)
                    polyhedra_kind['edges'] = polyhedra_kind['edges'] + list(edge)
                    polyhedra_kind['faces'] = polyhedra_kind['faces'] + list(face)
                    # print('edge: ', edge)
                    for e in edge:
                        # print(e)
                        center = (polyhedra_kind['vertices'][e[0]] + polyhedra_kind['vertices'][e[1]])/2.0
                        vec = polyhedra_kind['vertices'][e[0]] - polyhedra_kind['vertices'][e[1]]
                        length = np.linalg.norm(vec)
                        nvec = vec/length
                        # print(center, nvec, length)
                        polyhedra_kind['edge_cylinder']['lengths'].append(length/2.0)
                        polyhedra_kind['edge_cylinder']['centers'].append(center)
                        polyhedra_kind['edge_cylinder']['normals'].append(nvec)
            if len(polyhedra_kind['vertices']) > 0:
                polyhedra_kinds[spi] = polyhedra_kind
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
            if bondsetting[(spi, spj)][2] > 0:
                bondlists0 = np.append(bondlists0, bondlists1[speciesarray[bondlists1[:, 1]] == spj], axis = 0)
        bondlists1 = bondlistsji[speciesarray[bondlistsij[:, 0]] == spi]
        for spj in specieslist:
            if (spj, spi) not in bondsetting: continue
            if bondsetting[(spj, spi)][2] > 0:
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
    bondlists = build_bondlists(atoms, bondsetting)
    datas = build_polyhedralists(atoms, bondlists, bondsetting)
    