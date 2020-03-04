from ase.io import read, write
from blase.tools import write_blender, get_polyhedra_kinds, get_bondpairs
atoms = read('datas/anthraquinone.cif')
kwargs = {'show_unit_cell': 1, 
          'engine': 'BLENDER_WORKBENCH', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH', CYCLES
          'radii': 0.6,
          'bond_cutoff': 1.0,
          'display': True,
          'search_pbc': {'molecule_list':[['C', 'C'], ['C', 'O']]},
          'outfile': 'figs/test-search-molecule',
          }
write_blender(atoms, **kwargs)
