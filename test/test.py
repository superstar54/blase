import pytest
from blase.batoms import Batoms
import numpy as np


def test_build():
    from blase.build import surface, bulk
    pt = bulk('pt', 'Pt', cubic  =True)
    assert isinstance(pt, Batoms)
    pt111 = surface('pt111', pt, (1, 1, 1), 4, vacuum = 5)
    assert isinstance(pt111, Batoms)


def cavity():
    from blase.bio import read
    mof = read('docs/source/_static/datas/mof-5.cif')
    mof.bondsetting[('Zn', 'O')] = [0, 2.5, True, False]
    mof.bondsetting[('C', 'H')] = [0, 1.4, False, False]
    mof.boundary = 0.2
    mof.draw_cavity(9.0)
    mof.model_type = 2



def test_get_angles():
    from ase.build import molecule
    atoms = molecule('H2O')
    h2o = Batoms(atoms = atoms, label = 'h2o')
    angle = h2o.get_angle('H', 0, 'O', 0, 'H', 1)
    d = h2o.get_distances('H', 0, 'H', 1)


def test_extend():
    from ase.build import molecule
    co = molecule('CO')
    co = Batoms(label = 'co', atoms = co)
    co.translate([0, 0, 2])
    h2o = molecule('H2O')
    h2o = Batoms(label = 'h2o', atoms = h2o)
    h2o.extend(co)
