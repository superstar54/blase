import pytest
from blase.batoms import Batoms
from blase.butils import removeAll
import numpy as np


def test_batoms():
    """
    """
    removeAll()
    h2o = Batoms({'O': [[0, 0, 0.40]], 'H': [[0, -0.76, -0.2], [0, 0.76, -0.2]]})
    assert isinstance(h2o, Batoms)
    # properties
    h2o.cell = [3, 3, 3]
    h2o.pbc = True
    h2o.model_type = 1
    h2o.translate([0, 0, 2])
    #
    h2o_2 = h2o.copy('h2o_2')
    assert isinstance(h2o_2, Batoms)
    h2o_3 = Batoms(from_collection='h2o_2')
    #
    # delete
    del h2o_3['H'][[0]]
    assert len(h2o_3['H']) == 1
    #
    h2o_3.replace('O', 'C', [0])
    #
    h2o['O'][0] = [0, 0, 2.0]
    index = [0, 1]
    h2o['H'][index][:, 0] += 2

def test_canvas():
    from ase.build import fcc111
    from blase.batoms import Batoms
    from blase.butils import removeAll
    removeAll()
    atoms = fcc111('Pt', size = (4, 4, 4), vacuum=0)
    pt111 = Batoms(label = 'pt111', atoms = atoms)
    pt111.cell[2, 2] += 5
    pt111.render([1, 1, 1])

def test_batoms_animation():
    from ase.build import molecule
    removeAll()
    atoms = molecule('C2H6SO')
    images = []
    for i in range(20):
        temp = atoms.copy()
        temp.rotate(18*i, 'z')
        images.append(temp)
    c2h6so = Batoms(label = 'c2h6so', atoms = images)
    c2h6so.load_frames()
    # c2h6so.render(animation = True)


if __name__ == '__main__':
    test_batoms()
    test_batoms_animation()
    print('\n Batoms: All pass! \n')