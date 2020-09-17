from ase.build import molecule
from ase.io import read
from blase.tools import get_bondpairs, write_blender

atoms = read('datas/ATP.pdb')
atoms.positions[:, 2] -= min(atoms.positions[:, 2]) - 2
camera_loc = atoms[3].position + [50, 0, 20]

kwargs = {'show_unit_cell': 0, 
          'engine': 'BLENDER_WORKBENCH', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH'
          'camera_loc': camera_loc,  # distance from camera to front atom
          'camera_type': 'PERSP',  #  ['PERSP', 'ORTHO', 'PANO']
          'camera_lens': 100,  #
          'camera_target': atoms[3].position, #
          'ortho_scale': None, #
          'radii': 0.6, 
          'bond_cutoff': 1.0,
          'display': True,
          # 'functions': [['draw_plane', {'size': 1000, 'location': (0, 0, -1.0)}]],
          # 'outfile': 'figs/atp-cycles-1000',
          }
write_blender(atoms, **kwargs)