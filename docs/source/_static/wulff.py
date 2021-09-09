from ase.cluster import wulff_construction
from ase.visualize import view
from runblase import write_blender

surfaces = [(1, 1, 1), (1, 0, 0)]
energies = [1.28, 1.69]
atoms = wulff_construction('Au', surfaces, energies, 5000, 'fcc')
atoms.center(vacuum=2.0)
# view(atoms)
camera_loc =  atoms.get_center_of_mass() + [0, -300, 200]

batoms = {'atoms': atoms,
        'model_type': '0',
        'show_unit_cell': False,
        }
blase = {
          'engine': 'CYCLES', #'BLENDER_EEVEE', 'BLENDER_WORKBENCH', 'CYCLES'
          'camera_loc': camera_loc,  # distance from camera to front atom
          'camera_type': 'PERSP',  #  ['PERSP', 'ORTHO', 'PANO']
          'camera_lens': 100,  #
          'camera_target': atoms.get_center_of_mass(), #
          'functions': [['draw_plane', {'size': 1000, 'location': (0, 0, -1.0),  'color': [1.0, 1.0, 1.0, 1.0], 'material_style': 'mirror'}]],
          'output_image': 'figs/wulff',
  }
write_blender(batoms, blase, display = True)