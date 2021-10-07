from ase.io import read
from blase.butils import removeAll
from blase.batoms import Batoms

removeAll()
atoms = read('datas/tio2.cif')
tio2 = Batoms(label = 'tio2', atoms = atoms, color_style='VESTA')
tio2.boundary = 0.01
tio2.draw_cell()
for model_type in ['0', '1', '2', '3']:
    tio2.model_type = model_type
    tio2.render.run(direction = [0, 0, 1], engine = 'eevee', resolution_x = 500, output = 'batoms_model_type_%s.png'%model_type)
