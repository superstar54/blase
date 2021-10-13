from ase import atoms
from ase.io import read
from blase.batoms import Batoms
from blase.bdraw import draw_plane
from blase.butils import removeAll
removeAll()
cto = read('docs/source/_static/datas/catio3.cif')
cto = Batoms('cto', atoms = cto)
# ball
cto['Ca'].bsdf = {'Base Color': [0, 1, 0, 1.0], 'Metallic': 0.9, 'Specular': 1.0, 'Roughness': 0.01 }
cto['Ti'].bsdf = {'Base Color': [0.7, 0.8, 0, 1.0], 'Metallic': 0.2, 'Specular': 0.2, 'Roughness': 0.6 }
cto.bondsetting['Ti-O'].search = 0
cto.boundary = 0.01
cto.bondsetting.remove_bonds(['Ca', 'O'])
cto.polyhedrasetting['Ti'].color = [0.7, 0.4, 0.0, 1.0]
cto.render.set_world(color = [0.2, 0.2, 0.2, 1.0])
draw_plane(location = [0, 0, -cto['Ti'].size[0]], size = 200, color = (0.9, 0.9, 0.9, 1))
cto.render.light_type = 'POINT'
cto.render.lock_light_to_camera = False
cto.render.light_energy = 50000
cto.render.light_loc = [10, 4, 10]
cto.render.run([1, -0.3, 0.3], engine = 'cycles', output = 'gallery_catio3_ball.png')


# polyhedra
cto.boundary = 0.0
cto['O'].scale = 0.2
cto['Ti'].scale = 0.5
cto['Ca'].scale = 0.5
cto.boundary = 0.01
cto.bondsetting['Ti-O'].search = 1
cto.bondsetting['Ti-O'].style = '0'
cto.bondsetting['Ti-O'].width = 0.01
cto.polyhedrasetting['Ti'].color = [0.5, 0.0, 0.7, 0.2]
cto.draw_bonds()
cto.draw_polyhedras()
cto.render.set_world(color = [0.2, 0.2, 0.2, 1.0])
draw_plane(location = [0, 0, -2.5], size = 200, color = (0.9, 0.9, 0.9, 1))
cto.render.run([1, -0.3, 0.3], engine = 'cycles', output = 'gallery_catio3_polyhedra.png')




# bond
cto.boundary = 0.0
cto.model_type = 1
cto.repeat([2, 2, 2])
cto.cell = cto.cell[:]/2
del cto['O'][(cto['O'].positions[:, 2] > 4) | (cto['O'].positions[:, 0] > 4) | (cto['O'].positions[:, 1] > 4)]
del cto['Ca'][(cto['Ca'].positions[:, 2] > 4) | (cto['Ca'].positions[:, 0] > 4) | (cto['Ca'].positions[:, 1] > 4)]
cto.bondsetting['Ti-O'].search = 0
cto.bondsetting['Ti-O'].color1 = [0.1, 0.1, 0.1, 1.0]
cto.model_type = 1
cto['O'].scale = 0.2
draw_plane(location = [0, 0, -cto['Ti'].size[0]], size = 200, color = (0.9, 0.9, 0.9, 1))
cto.render.run([1, -0.3, 0.3], engine = 'eevee',  output = 'gallery_catio3_bond.png')

