from ase.build import molecule
from blase.tools import get_bondpairs, write_blender
atoms = molecule('H2O')
atoms.center(vacuum=0.2)
atoms.rotate(90, 'y')

kwargs = {'show_unit_cell': 1, 
          'radii': 0.6, 
          'bond_cutoff': 1.0,
          # 'display': True,
          'outfile': 'figs/h2o',
          }


from blase.bio import Blase

bobj = Blase(atoms, **kwargs)
bobj.draw()
# bobj.render()
# bobj.export('h2o.xyz')
print('\n Finished!')