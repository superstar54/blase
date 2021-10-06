from ase.build import bulk
from blase.batoms import Batoms
au = bulk('Au', 'fcc', cubic=True)
au = Batoms(label = 'au', atoms = au)
au.draw_cell()
for light in ["Default", "basic.sl", "outdoor.sl", "paint.sl", "rim.sl", "studio.sl",]:
    au.render.run(direction = [1, -0.3, 0.1], engine = 'workbench', studiolight = light, sresolution_x = 200, output = 'studiolight_%s.png'%light)
