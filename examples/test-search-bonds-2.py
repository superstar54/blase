from ase.io import read, write
from ase.visualize import view
from blase.tools import write_blender, get_polyhedra_kinds, get_bondpairs
import numpy as np
from pprint import pprint
from ase.data import covalent_radii


atoms = read('datas/perovskite.cif')
atoms = atoms*[2, 2, 2]
# view(atoms)

kwargs = {'show_unit_cell': 1, 
          'engine': 'BLENDER_WORKBENCH', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH', CYCLES
          'radii': 0.6,
          'bond_cutoff': 1.0,
          'display': True,
          'search_pbc': {'bonds_dict': {'I': [['Pb'], -1]}, 'molecule_list':[['C', 'N']]},
          'polyhedra_dict': {'Pb': ['I']},
          'outfile': 'figs/test-search-bonds-2',
          }
write_blender(atoms, **kwargs)
