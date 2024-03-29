import pytest
from blase.butils import removeAll
from blase.batoms import Batoms
from blase.bio import read
import numpy as np



def test_boundary():
    from blase.butils import removeAll
    from blase.batoms import Batoms
    removeAll()
    a = 3.96
    positions = [[0, 0, 0], [a/2, a/2, 0], [a/2, 0, a/2], [0, a/2, a/2]]
    pt = Batoms('pt', {'Pt': positions}, pbc = True, cell = (a, a, a))
    pt.boundary = 0.01
    atoms = pt.get_atoms_with_boundary()
    assert len(atoms) == 14


def test_bond_search():
from blase.butils import removeAll
from blase.batoms import Batoms
from blase.bio import read
removeAll()
mol = read('docs/source/_static/datas/c2h6so.xyz')
mol.cell = [4, 4, 4]
mol.pbc = True
mol.boundary = 0.01
mol.model_type = 1

def test_wrap():
    from blase.butils import removeAll
    from blase.batoms import Batoms
    from blase.bio import read
    removeAll()
    tio2 = read('../docs/source/_static/datas/tio2.cif')
    tio2.boundary = 0.01
    tio2.model_type = 2
    atoms = tio2.get_atoms_with_boundary()
    assert len(atoms) == 55



if __name__ == '__main__':
    test_boundary()
    test_bond_search()
    print('\n Bondsetting: All pass! \n')