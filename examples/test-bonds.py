from ase.io import read, write
from ase.visualize import view
from runblase import write_blender
import numpy as np
from pprint import pprint
from ase.data import covalent_radii


atoms = read('datas/tio2.cif')
atoms = atoms*[2, 2, 2]
atoms.pbc = False
# view(atoms)
# pprint(bond_list)
kwargs = {
    'name': 'tio2',
        #   'show_unit_cell': 1, 
          'engine': 'BLENDER_WORKBENCH', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH', CYCLES
        #   'radii': 0.4,
        #   'bond_cutoff': 1.0,
          'display': True,
          'outfile': 'figs/test-search-bonds',
          }
write_blender(atoms, **kwargs)
