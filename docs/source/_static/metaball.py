from ase.cluster import wulff_construction
from ase.build import molecule
from ase.visualize import view
from blase.tools import get_bondpairs, write_blender

surfaces = [(1, 1, 1), (1, 0, 0)]
energies = [1.28, 1.69]
atoms = wulff_construction('Au', surfaces, energies, 10, 'fcc')
atoms.cell = [10, 10, 10]
kwargs = {
		  'show_unit_cell': 1, 
		  'balltypes': {'Au': 'meta'},
          'radii': 2.0,
          'outfile': 'figs/test-metaball',
          }
write_blender(atoms, **kwargs)

