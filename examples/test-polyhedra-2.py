from ase.io import read, write
from blase.tools import write_blender

atoms = read('datas/tio2.cif')
atoms = atoms*[2, 2, 2]
kwargs = {'show_unit_cell': 1, 
          'engine': 'BLENDER_WORKBENCH',
          'radii': 0.6,
          'bond_cutoff': 1.0,
          'search_pbc_atoms': {'bonds_dict': {'O': [['Ti'], -1]}},
          'polyhedra_dict': {'Ti': ['O']},
          'display': True,
          # 'outfile': 'figs/test-search-bonds',
          }
write_blender(atoms, **kwargs)