from ase.build import molecule
from blase.tools import get_bondpairs, write_blender

atoms = molecule('CH3OH')

kwargs = {'show_unit_cell': 0, 
          'engine': 'BLENDER_WORKBENCH', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH'
          'radii': 0.6, 
          'bond_cutoff': 0.8,
          'display': True,
          # 'functions': [['highlight_atoms', {'atomlist': [0]}]],
          # 'outfile': 'ch3oh',
          'debug': True,
          }
write_blender(atoms, **kwargs)
