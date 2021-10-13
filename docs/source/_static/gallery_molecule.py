# Ball
from ase.build import molecule
from blase.batoms import Batoms
from blase.bdraw import draw_plane
from blase.butils import removeAll
removeAll()
h2o = molecule('H2O')
h2o = Batoms('h2o', atoms = h2o)
h2o['H'].color = [0.8, 1.0, 0.0, 1.0]
draw_plane(location = [0, 0, h2o['H'][0][2] - h2o['H'].size[0]], size = 400, color = (0.1, 0.1, 0.1, 1))
h2o.render.light_type = 'POINT'
h2o.render.lock_light_to_camera = False
h2o.render.light_energy = 50000
h2o.render.light_loc = [2, 0, 10]
h2o.render.camera_loc = [4, 0, 0.3]
h2o.render.run(engine = 'cycles', ortho_scale = 2.5, ratio = 0.8, output = 'gallery_h2o_ball.png')


# Bond
from ase.build import molecule
from blase.batoms import Batoms
from blase.bdraw import draw_plane
from blase.butils import removeAll
removeAll()
h2o = molecule('H2O')
h2o = Batoms('h2o', atoms = h2o)
h2o['O'].scale = 0.6
h2o['H'].scale = 0.6
h2o['H'].color = [0.8, 1.0, 0.0, 1.0]
h2o.bondsetting['O-H'].width = 0.1
h2o.bondsetting['O-H'].style = '0'
h2o.bondsetting['O-H'].color1 = [0.8, 0.8, 0.8, 1.0]
h2o.draw_bonds()
draw_plane(location = [0, 0, h2o['H'][0][2] - h2o['H'].size[0]], size = 400, color = (0.1, 0.1, 0.1, 1))
h2o.render.light_type = 'POINT'
h2o.render.lock_light_to_camera = False
h2o.render.light_energy = 50000
h2o.render.light_loc = [2, 0, 10]
h2o.render.camera_loc = [4, 0, 0.3]
h2o.render.run(engine = 'cycles', ortho_scale = 2.5, ratio = 0.8, output = 'gallery_h2o_bond.png')
