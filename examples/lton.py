from ase.io import read, write
from ase.visualize import view
from blase.tools import write_blender, get_polyhedra_kinds, get_bondpairs
import numpy as np
from pprint import pprint
from ase.data import covalent_radii


atoms = read('datas/latino2-trans-p-001-plp-relax.in')
# atoms.pbc = [False, False, False]
# atoms = atoms*[2, 2, 2]
# kind_props = {
# 'Ti': {'radius': 0.6, 'color': [0/255.0, 191/255.0, 56/255.0]},
# 'O': {'radius': 0.6, }
# }

bond_list = get_bondpairs(atoms, cutoff=1.0, rmbonds=[['Ti', 'Ti'], ['La', 'O'], ['La', 'N'], ['La', 'Ti'], ['La', 'La']])

# bond_list = get_bondpairs(atoms, rmbonds = [['Ti', 'Ti']])
kwargs = {'show_unit_cell': 1, 
          'engine': 'BLENDER_WORKBENCH', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH', CYCLES
          'radii': 0.3,
          # 'bond_cutoff': 1.0,
          'bond_list': bond_list,
          # 'kind_props': kind_props,
          'display': True,
          'polyhedra_dict': {'Ti': ['O', 'N']},
          # 'search_pbc_atoms': {'bonds_dict':{'I': [['Pb'], -1]}},
          # 'outfile': 'figs/test-polyhedra-lton',
          }
write_blender(atoms, **kwargs)