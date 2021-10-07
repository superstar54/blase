from ase.build import bulk
from blase.batoms import Batoms
au = bulk('Au', 'fcc', cubic=True)
au = Batoms(label = 'au', atoms = au)
au.draw_cell()
for engine in ['workbench', 'eevee', 'cycles']:
    au.render.run(direction = [1, -0.3, 0.1], engine = engine, resolution_x = 200, output = 'render_%s.png'%engine)
