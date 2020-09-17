from ase.build import molecule
from ase.visualize import view
from blase.tools import get_bondpairs, write_blender

atoms = molecule('C2H6SO')
atoms.center(vacuum  = 1.0)
camera_loc = atoms[3].position + [30, 0, 10]

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
          # 'functions': [['draw_plane', {'size': 100, 'location': (0, 0, -0.0), 'color': (0.1, 0.1, 0.1, 1.0), 'material_style': 'mirror'}]],
          # 'outfile': 'figs/c2h6so-cycles',
          }
write_blender(atoms, **kwargs)