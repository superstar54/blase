from ase.build import molecule
from ase.io import read
from ase.visualize import view
from blase.tools import get_bondpairs, write_blender

atoms = read('datas/sco-1-6ml-oh.pwo')
# print(atoms)
atoms.write('datas/sco-1-6ml-oh.in')
atoms = read('datas/sco-1-6ml-oh.in')
atoms = atoms*[2, 2, 1]
print(atoms)

kwargs = {'show_unit_cell': 1, 
          'engine': 'BLENDER_WORKBENCH', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH'
          'radii': 0.6, 
          'bond_cutoff':  0.8,
          'display': True,
          'polyhedra_dict': {'Cr': ['O', 'N']},
          # 'debug': True,
          # 'outfile': 'figs/c2h6so-cycles',
          }
write_blender(atoms, **kwargs)