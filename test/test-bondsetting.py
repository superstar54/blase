import pytest
from blase.butils import removeAll
from blase.batoms import Batoms
from blase.bio import read
import numpy as np


def test_replace():
    from ase.build import molecule
    from blase.batoms import Batoms
    from blase.butils import removeAll
    removeAll()
    co = Batoms('co', atoms = molecule('CO'))
    co.cell = [2, 2, 3]
    co.pbc = True
    co.repeat([2, 2, 2])
    co.replace('C', 'C_1', [5])
    co.model_type = 1
    

def test_boundary():
    a = 3.96
    positions = [[0, 0, 0], [a/2, a/2, 0], [a/2, 0, a/2], [0, a/2, a/2]]
    pt = Batoms({'Pt': positions}, pbc = True, cell = (a, a, a))
    pt.boundary = 0.01


def test_polyhedra():
    tio2 = read('docs/source/_static/datas/tio2.cif')
    tio2.model_type = 2

def test_search_bond():
    pk = read('docs/source/_static/datas/perovskite.cif')
    pk.model_type = 2

if __name__ == '__main__':
    test_replace()
    test_boundary()
    test_polyhedra()
    print('\n Bondsetting: All pass! \n')