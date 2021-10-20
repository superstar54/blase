from ase.io import read
from blase.batoms import Batoms
from blase.bdraw import draw_plane
from blase.butils import removeAll
removeAll()
atoms = read('docs/source/_static/datas/mof-5.cif')
mof = Batoms(label = 'mof-5', atoms = atoms)
mof['H'].color = [0.6, 0, 1.0, 1.0]
mof['C'].color = [0.0, 0.6, 0.1, 1.0]
mof.polyhedrasetting['Zn'].color = [0.1, 0.4, 0.7, 1.0]
mof.model_type = 2
mof.draw_cavity_sphere(9.0, boundary = [[0.2, 0.8], [0.2, 0.8], [0.2, 0.8]])
mof.draw_cell()
mof.render.set_world(color = [0.2, 0.2, 0.2, 1.0])
draw_plane(location = [0, 0, 0], size = 200, color = (0.9, 0.9, 0.9, 1))
mof.render.light_type = 'POINT'
mof.render.lock_light_to_camera = False
mof.render.light_energy = 500000
mof.render.light_loc = [40, 4, 40]
mof.render.run([1, -0.3, 0.3], engine = 'cycles', output = 'gallery_cavity.png')