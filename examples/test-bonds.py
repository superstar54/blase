from ase.io import read, write
from ase.visualize import view
from runblase import write_blender
import numpy as np
from pprint import pprint
from ase.data import covalent_radii


atoms = read('datas/tio2.cif')
atoms = atoms*[2, 1, 1]
atoms.pbc = False
atoms.info['species'] = atoms.get_chemical_symbols()
atoms.info['species'][0] = 'Ti_1'
# view(atoms)
# pprint(bond_list)

batoms = {
        'model_type': '2',
        # 'show_unit_cell': False,
        'remove_bonds': {'Ti_1':['O']},
        'boundary': [0.00, 0.00, 0.00],
        'polyhedra_dict': {'Ti': ['O', 'N']},
        'color': 'VESTA',
        }
blase = {
          # 'engine': 'BLENDER_WORKBENCH', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH', CYCLES
          'output_image': 'figs/test-search-bonds-1',
          'display': True,
  }

write_blender(atoms, batoms, blase)
