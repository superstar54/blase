from ase.io import read, write
from runblase import write_blender

atoms = read('datas/ATP.pdb')
atoms.positions[:, 2] -= min(atoms.positions[:, 2]) - 2
camera_loc = atoms[3].position + [50, 0, 20]
batoms = {'atoms': atoms, 'model_type': '1'}
blase = {
          'engine': 'CYCLES', #'BLENDER_EEVEE', 'BLENDER_WORKBENCH', 'CYCLES'
          'camera_loc': camera_loc,  # distance from camera to front atom
          'camera_type': 'PERSP',  #  ['PERSP', 'ORTHO', 'PANO']
          'camera_lens': 100,  #
          'camera_target': atoms[3].position, #
          'functions': [['draw_plane', {'size': 1000, 'location': (0, 0, -1.0),  'color': [1.0, 1.0, 1.0, 1.0]}]],
          'output_image': 'figs/atp',
  }
write_blender(batoms, blase)

