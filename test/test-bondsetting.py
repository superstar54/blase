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
    pt = Batoms('pt', {'Pt': positions}, pbc = True, cell = (a, a, a))
    pt.boundary = 0.01


def test_polyhedra():
    tio2 = read('../docs/source/_static/datas/tio2.cif')
    tio2.model_type = 2

def test_search_bond():
    pk = read('../docs/source/_static/datas/perovskite.cif')
    pk.model_type = 2

def test_hydrogen_bond():
    from ase.build import molecule
    from blase.batoms import Batoms
    h2o = molecule('H2O')
    h2o2 = molecule('H2O')
    h2o2.rotate(90, 'x')
    h2o2.translate([0, 0, 3])
    h2o = h2o + h2o2
    h2o = Batoms(label = 'h2o', atoms = h2o)
    h2o.bondsetting['H-O'].min = 2.0
    h2o.bondsetting['H-O'].max = 3.0
    h2o.bondsetting['H-O'].bondlinewidth = 0.01
    h2o.model_type = 1
    h2o.render.run([1, 0 ,0], engine = 'eevee')



if __name__ == '__main__':
    test_replace()
    test_boundary()
    test_polyhedra()
    test_hydrogen_bond()
    print('\n Bondsetting: All pass! \n')