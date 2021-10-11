from ase.build import molecule
from blase.butils import removeAll
from blase.batoms import Batoms

removeAll()
co = Batoms(label = 'co', atoms = molecule('CO'))
co.bondsetting['C-O'].style = '0'
co.model_type = 1
co.render.run(direction = [1, 0, 0], engine = 'eevee', resolution_x = 300, output = 'bondsetting_style_0.png')
#
co.bondsetting['C-O'].style = '1'
co.model_type = 1
co.render.run(direction = [1, 0, 0], engine = 'eevee', resolution_x = 300, output = 'bondsetting_style_1.png')
#
co.bondsetting['C-O'].style = '2'
co.bondsetting['C-O'].width = 0.01
co.model_type = 1
co.render.run(direction = [1, 0, 0], engine = 'eevee', resolution_x = 300, output = 'bondsetting_style_2.png')
#
co.bondsetting['C-O'].style = '3'
co.bondsetting['C-O'].width = 0.01
co.model_type = 1
co.render.run(direction = [1, 0, 0], engine = 'eevee', resolution_x = 300, output = 'bondsetting_style_3.png')
