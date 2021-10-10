"""
"""
import collections
import bpy
from blase.butils import object_mode
import numpy as np
from time import time
from pprint import pprint


class Setting():
    """
    Setting object

    The Setting object store the other information for batoms.

    Parameters:

    label: str
        The label define the batoms object that a Setting belong to.
    """
    def __init__(self, label) -> None:
        self.label = label
        self.name = 'base'
    def get_data(self):
        data = {}
        for b in self.collection:
            data[b.name] = b
        return data
    @property
    def collection(self):
        return self.get_collection()
    def get_collection(self):
        collection = getattr(bpy.data.collections[self.label], self.name)
        return collection
    @property
    def species(self):
        return self.get_species()
    def get_species(self) -> dict:
        """
        read species from collection.
        """
        species = {}
        coll_atom = bpy.data.collections['%s_atom'%self.label]
        for ba in coll_atom.objects:
            species[ba.blasebatom.species] = {'color': ba.children[0].data.materials[0].diffuse_color,'radius': ba.blasebatom.radius}
        return species
    @property
    def data(self):
        return self.get_data()
    def __getitem__(self, index):
        item = self.find(index)
        if item is None:
            raise Exception('%s not in %s setting'%(index, self.name))
        return item
    def __setitem__(self, index, value):
        """
        Add bondpair one by one
        """
        collection = self.collection
        item = self.find(index)
        
    def copy(self, label):
        object_mode()
        bondsetting = self.__class__(label)
        for key, b in self.data.items():
            bondsetting[key] = b.as_list()
        return bondsetting
    def __add__(self, other):
        self += other
        return self
    def __iadd__(self, other):
        self.extend(other)
        return self
    def extend(self, other):
        for key, value in other.data.items():
            self[key] = value.as_list()
        return self
    def __iter__(self):
        item = self.collection
        for i in range(len(item)):
            yield item[i]
    def __len__(self):
        return len(self.collection)
    def find(self, name):
        i = self.collection.find(name)
        if i == -1:
            return None
        else:
            return self.collection[i]


class Bondsetting(Setting):
    """
    Bondsetting object

    The Bondsetting object store the bondpair information.

    Parameters:

    label: str
        The label define the batoms object that a Bondsetting belong to.
    cutoff: float
        Cutoff used to calculate the maxmium bondlength for bond pairs.
    """
    def __init__(self, label, cutoff = 1.3) -> None:
        Setting.__init__(self, label)
        self.label = label
        self.name = 'blasebond'
        self.cutoff = cutoff
        self.set_default(self.species, cutoff)
    def __setitem__(self, index, value):
        """
        Add bondpair one by one
        """
        bond = self.find(index)
        if bond is None:
            bond = self.collection.add()
        bond.symbol1 = index.split('-')[0]
        bond.symbol2 = index.split('-')[1]
        bond.name = index
        bond.min = value[0]
        bond.max = value[1]
        bond.search = value[2]
        bond.polyhedra = value[3]
        if len(value) == 8:
            bond.color1 = value[4]
            bond.color2 = value[5]
            bond.bondlinewidth = value[6]
            bond.style = value[7]
    def set_default(self, species, cutoff = 1.3):
        """
        """
        bondtable = get_bondtable(species, cutoff=cutoff)
        for key, value in bondtable.items():
            self['%s-%s'%(key[0], key[1])] = value
    def extend(self, other):
        for key, value in other.data.items():
            self[key] = value.as_list()
        # new
        species1 = set(self.species)
        species2 = set(other.species)
        nspecies1 = species1 - species2
        nspecies2 = species2 - species1
        for sp1 in nspecies1:
            for sp2 in nspecies2:
                species = {sp1: self.species[sp1], sp2: other.species[sp2]}
                self.set_default(species, self.cutoff)
    def add_bonds(self, bondpair):
        for key in bondpair:
            species = {sp: self.species[sp] for sp in key}
            self.set_default(species)
    def remove_bonds(self, bondpair):
        for key in bondpair:
            name = '%s-%s'%(key[0], key[1])
            i = self.collection.find(name)
            if i != -1:
                self.collection.remove(i)
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Bondpair      min     max   Search_bond    Polyhedra \n'
        for b in self.collection:
            s += '{0:10s} {1:4.3f}   {2:4.3f}      {3:10s}   {4:10s} \n'.format(\
                b.name, b.min, b.max, str(b.search), str(b.polyhedra))
        s += '-'*60 + '\n'
        return s


def get_bondtable(speciesdict, cutoff = 1.3):
    """
    """
    from blase.data import default_bonds
    bondtable = {}
    for species1 in speciesdict:
        element1 = species1.split('_')[0]
        color1 = speciesdict[species1]['color']
        radius1 = cutoff * speciesdict[species1]['radius']
        for species2 in speciesdict:
            element2 = species2.split('_')[0]
            pair = (element1, element2)
            if pair not in default_bonds: continue
            color2 = speciesdict[species2]['color']
            radius2 = cutoff * speciesdict[species2]['radius']
            bondmax = radius1 + radius2
            bondtable[(species1, species2)] = [0.5, bondmax, default_bonds[pair][0], default_bonds[pair][1], color1, color2, 0.10, '0']
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
    for b in bondsetting:
        spi = b.symbol1
        spj = b.symbol2
        eli = spi.split('_')[0]
        elj = spj.split('_')[0]
        key = (eli, elj)
        cutoff_min[key] = b.min
        cutoff[key] = b.max
    #
    tstart = time()
    nli, nlj, nlS = neighbor_list('ijS', atoms, cutoff, cutoff_min, self_interaction=False)
    bondlists = np.append(np.array([nli, nlj], dtype=int).T, np.array(nlS, dtype=int), axis = 1)
    # print('build_bondlists: {0:10.2f} s'.format(time() - tstart))
    return bondlists


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
            name = '%s-%s'%(spi, spj)
            if not bondsetting.find(name): continue
            if bondsetting[name].search > 0:
                bondlists0 = np.append(bondlists0, bondlists1[speciesarray[bondlists1[:, 1]] == spj], axis = 0)
        bondlists1 = bondlistsji[speciesarray[bondlistsji[:, 0]] == spi]
        for spj in specieslist:
            name = '%s-%s'%(spj, spi)
            if not bondsetting.find(name): continue
            if bondsetting[name].search > 0:
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
    