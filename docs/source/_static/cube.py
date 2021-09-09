from ase.io.cube import read_cube_data
from runblase import write_blender

data, atoms = read_cube_data('datas/h2o-homo.cube')
batoms = {'atoms': atoms,
          'model_type': '1',
          'show_unit_cell': True, 
          'isosurface': [data, -0.002, 0.002],
          }
blase = {
          'engine': 'BLENDER_WORKBENCH', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH'
          'bbox': [[0, 10], [0, 10], [0, 10]], # set range for the box
          'output_image': 'figs/h2o-homo-cube',
          }
write_blender(batoms, blase) 