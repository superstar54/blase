from ase.io import read
from blase.batoms import Batoms

atoms = read('docs/source/_static/datas/tio2.cif')
tio2 = Batoms(label = 'tio2', atoms = atoms, color_style='VESTA')
tio2.bondsetting[('Ti', 'O')] = [0.5, 2.5, True, True]
tio2.boundary = 0.01
tio2.draw_cell()
for model_type in ['0', '1', '2', '3']:
    tio2.model_type = model_type
    tio2.render.run(direction = [0, 0, 1], engine = 'eevee', resolution_x = 500, output = 'docs/source/_static/batoms_model_type_%s.png'%model_type)
