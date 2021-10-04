import pytest
from blase.render import Render
from blase.batoms import Batoms
from blase.butils import removeAll
import numpy as np


def test_render():
    """
    """
    removeAll()
    r = Render('blase')
    assert isinstance(r, Render)
    #
    from ase.build import fcc111
    removeAll()
    atoms = fcc111('Au', size = (4, 4, 4), vacuum=0)
    au111 = Batoms(label = 'au111', atoms = atoms)
    au111.cell[2, 2] += 5
    au111.render([1, 1, 1], output_image='au111')
    au111.render([1, 1, 1], engine = 'eevee', output_image='au111-eevee')
    au111.render([1, 1, 1], engine = 'cycles', output_image='au111-cycles')
    print(au111.rendersetting)

test_render()