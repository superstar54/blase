from ase.io.cube import read_cube_data
from blase.tools import write_blender

data, atoms = read_cube_data('test.cube')
kwargs = {'show_unit_cell': 1, 
          'engine': 'BLENDER_EEVEE', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH'
          'radii': 0.4, 
          # 'bond_cutoff': 1.0,
          'isosurface': [data, -0.002, 0.002],
          'rotations': [[90, '-x'], [30, 'y']],
          'display': True,
          'outfile': 'figs/testcube'}

write_blender(atoms, **kwargs)