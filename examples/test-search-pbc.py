from ase.io import read, write
from ase.build import bulk
from ase.visualize import view
from blase.tools import write_blender, get_polyhedra_kinds, get_bondpairs
import numpy as np
from pprint import pprint
from ase.data import covalent_radii

atoms = bulk('Pt')
kwargs = {'show_unit_cell': 1, 
          # 'engine': 'BLENDER_WORKBENCH', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH', CYCLES
          'radii': 0.6,
          'display': True,
          'search_pbc_atoms': {'bonds_dict': {}},
          'outfile': 'figs/test-search-pbc',
          }
write_blender(atoms, **kwargs)
