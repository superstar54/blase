from ase.io import read
from blaseio import write_blender

atoms = read('datas/po.xyz')
batoms = {'atoms': atoms, 'model_type': '1'}
camera_loc = atoms[0].position + [20, 5, 20]
blase = {
          'engine': 'CYCLES', #'BLENDER_EEVEE', 'BLENDER_WORKBENCH', 'CYCLES'
          'camera_loc': camera_loc,
          'camera_type': 'PERSP', 
          'camera_lens': 100,  #
          'camera_target': atoms.get_center_of_mass(), #
          'functions': [['draw_plane', {'size': 100, 'location': (0, 0, -1.0),  'color': [1.0, 1.0, 1.0, 1.0]}]],
          'output_image': 'figs/po',
        #   'run_render': False,
  }
for model_type in ['0', '1', '2', '3']:
    batoms['model_type'] = model_type
    blase['output_image'] = 'po-%s.png'%model_type
    write_blender(batoms, blase)
