from ase.io.cube import read_cube_data
from blase.batoms import Batoms
volume, atoms = read_cube_data('docs/source/_static/datas/h2o-homo.cube')
h2o = Batoms('h2o', atoms = atoms, volume = volume, draw = False)
h2o.isosurfacesetting[1] = [0.002, [1, 1, 0, 0.5]]
h2o.isosurfacesetting[2] = [-0.002, [0, 0, 0.8, 0.5]]
h2o.draw_isosurface()
