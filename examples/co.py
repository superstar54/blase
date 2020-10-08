from ase.build import molecule
from ase.atoms import Atoms
from ase.visualize import view
from blase.tools import get_bondpairs, write_blender
import time

atoms = Atoms('CO', positions = [[0, 0, 1e-6,], [1.5, 1e-6, 0]])
kwargs = {'show_unit_cell': 0, 
          'engine': 'BLENDER_WORKBENCH', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH'
          'radii': 0.6, 
          'bond_cutoff':  1.0,
          # 'display': True,
          'run_render': False,
          'outfile': 'figs/h2o',
          }
tstart = time.time()
write_blender(atoms, **kwargs)
print('run: {0:10.2f} s'.format(time.time() - tstart))

