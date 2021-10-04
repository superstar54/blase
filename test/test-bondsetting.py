import pytest
from blase.batoms import Batoms
from blase.bio import read
import numpy as np


def test_boundary():
    a = 3.96
    positions = [[0, 0, 0], [a/2, a/2, 0], [a/2, 0, a/2], [0, a/2, a/2]]
    pt = Batoms({'Pt': positions}, pbc = True, cell = (a, a, a))
    pt.boundary = 0.01


def test_polyhedra():
    tio2 = read('docs/source/_static/datas/tio2.cif')
    tio2.bondsetting[('Ti', 'O')] = [0.5, 2.5, True, False]
    tio2.model_type = 2

if __name__ == '__main__':
    test_boundary()
    test_polyhedra()
    print('\n Bondsetting: All pass! \n')