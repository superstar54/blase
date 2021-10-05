from ase.build import bulk
from blase.batoms import Batoms
from ase.io import read
au = bulk('Au', cubic = True)
au = Batoms(label = 'au', atoms = au, segments = [6, 6])
au.repeat([50, 50, 100])
au.render.run(direction = [1, 0, 0])
