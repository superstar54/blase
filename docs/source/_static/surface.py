from ase.build import fcc111
from blase.batoms import Batoms
import numpy as np

atoms = fcc111('Au', size = (8, 8, 4), vacuum=0)
au111 = Batoms(label = 'au111', atoms = atoms)
au111.cell[2, 2] += 10
au111.draw_cell()
au111.render.run(direction = [0, 0, 1], engine = 'eevee', camera_type = 'ORTHO', output = 'docs/source/_static/gallery_top_view.png')
au111.render.run(direction = [1, 0, 0], engine = 'eevee', camera_type = 'ORTHO', output = 'docs/source/_static/gallery_side_view.png')
au111.render.run(direction = [1, -0.3, 0.1],
                engine = 'eevee', 
                camera_type = 'PERSP',
                camera_loc = au111['Au'][227] + np.array([-50, 0, 10]),
                ratio = 1,
                output = 'docs/source/_static/gallery_persp_view.png')