from ase.io import read, write
from ase.visualize import view
from blase.tools import write_blender, get_polyhedra_kinds, get_bondpairs, search_pbc, search_molecule
import numpy as np
from pprint import pprint
from ase.data import covalent_radii


atoms = read('perovskite.cif')
# view(atoms)
newatoms = search_pbc(atoms, remove_dict = {'I': ['Pb']})
newatoms = search_molecule(newatoms, search_list = ['C', 'N'])
view(newatoms)


# bond_list = get_bondpairs(atoms, rmbonds = [['Ti', 'Ti']])
kwargs = {'show_unit_cell': 1, 
          'engine': 'BLENDER_WORKBENCH', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH', CYCLES
          'radii': 0.6,
          'bond_cutoff': 1.0,
          # 'bond_list': bond_list,
          # 'kind_props': kind_props,
          'display': True,
          # 'world': True,
          'polyhedra_dict': {'Pb': ['I']},
          'outfile': 'figs/test-search-bonds-2',
          }
# write_blender(newatoms, **kwargs)

'''
from blase.bio import Blase

bobj = Blase(atoms, **kwargs)
bobj.draw_cell()
bobj.draw_atoms()
bobj.draw_bonds()
bobj.draw_polyhedras()
# bobj.render()
# bobj.export('h2o.xyz')
print('\n Finished!')
'''