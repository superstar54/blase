from ase.build import molecule
from blase.tools import get_bondpairs, write_blender

atoms = molecule('CH3OH')

kwargs = {'show_unit_cell': 0, 
          'radii': 0.6, 
          'bond_cutoff': 1.0,
          'display': True,
          # 'functions': [['highlight_atoms', {'atomlist': [0]}]],
          # 'outfile': 'ch3oh',
          }
write_blender(atoms, **kwargs)