from ase.build import molecule
from blase.tools import get_bondpairs, write_blender
atoms = molecule('H2O')
atoms.center(vacuum=1.0)
atoms.rotate(90, 'y')
# atoms = atoms*[5, 5, 5]

kwargs = {'show_unit_cell': 0, 
          'radii': 0.6, 
          'bond_cutoff': 1.0,
          # 'display': True,
          'outfile': 'figs/h2o',
          }
write_blender(atoms, **kwargs)

