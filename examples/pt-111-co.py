from ase.build import fcc111, add_adsorbate
from ase.atoms import Atoms
from ase.build import bulk
from ase.visualize import view
from blase.tools import get_bondpairs, write_blender

bulk = bulk('Pt')
bulk.write('pt.in')
adsorbate = Atoms('CO')
adsorbate[1].z = 1.1
atoms = fcc111('Pt', (2, 2, 3), a=3.96, vacuum=7.0)
add_adsorbate(atoms, adsorbate, 1.8, 'ontop')
# atoms = atoms*[5, 5, 1]
# view(atoms)

kwargs = {'show_unit_cell': 0, 
          'radii': 0.6,
          'bond_cutoff': 1.0,
          'display': True,
          'make_real': True,
          'outfile': 'figs/pt-111-co',
          }
write_blender(atoms, **kwargs)

