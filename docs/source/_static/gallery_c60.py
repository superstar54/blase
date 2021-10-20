from ase.build import molecule
from blase.batoms import Batoms
from blase.bdraw import draw_plane
from blase.butils import removeAll
import numpy as np
removeAll()
c60 = molecule('C60')
pos = c60.get_center_of_mass()
c60 = Batoms(label = 'c60', atoms = c60)
c60['C'].color = [0.1, 0.1, 0.1, 1.0]
c60.bondsetting['C-C'].color1 = [0.2, 0.8, 0.1, 1.0]
c60.bondsetting['C-C'].type = '0'
# we add a ghost site (species ``X``) at the center of a cavity
x = Batoms('x', {'X':[pos]})
c60 = c60 + x
# add bond `X-C`, and set polyhedra to True
c60.bondsetting['X-C'] = [0, 10, 2, True]
c60.polyhedrasetting['X'].color = [0.4, 0.4, 0, 1.0]
c60.model_type = 1
# # draw polyhedral model manually and not show the edge
c60.draw_polyhedras(show_edge = False)
c60.render.set_world(color = [0.2, 0.2, 0.2, 1.0])
draw_plane(location = [0, 0, min(c60.positions[:, 2]) - c60['C'].size[0]], 
            size = 200, color = (0.9, 0.9, 0.9, 1))
c60.render.light_type = 'POINT'
c60.render.lock_light_to_camera = False
c60.render.light_energy = 500000
c60.render.light_loc = [20, 4, 20]
c60.render.run([1, -0.3, 0.3], canvas = np.array([[-6, -6, -6], [6, 6, 6]]), 
        engine = 'cycles', output = 'gallery_c60.png')