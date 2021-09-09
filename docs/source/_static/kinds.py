from ase.build import fcc111
from runblase import write_blender
#============================================================
atoms = fcc111('Pt', (7, 7, 3), vacuum=0.0)
kind_props = {
'Pt_0': {'color': [208/255.0, 208/255.0, 224/255.0]},
'Pt_1': {'color': [225/255.0, 128/255.0, 0/255.0]},
'Pt_2': {'color': [0/255.0, 191/255.0, 56/255.0]},
}
atoms.info['species'] = []
for i in range(len(atoms)):
     ind = int((atoms[i].x/5))
     kind = atoms[i].symbol + '_{0}'.format(ind)
     atoms.info['species'].append(kind)

camera_loc = atoms.get_center_of_mass() + [0, -60, 30]

batoms = {'atoms': atoms,
          'kind_props': kind_props,
        'model_type': '0',
        'color': 'VESTA',
        }
blase = {
          'engine': 'CYCLES',
          'camera_loc': camera_loc,  # distance from camera to front atom
          'camera_type': 'PERSP',  #  ['PERSP', 'ORTHO', 'PANO']
          'camera_lens': 100,  #
          'camera_target': atoms.get_center_of_mass(), #
          'output_image': 'figs/kinds',
  }
write_blender(batoms, blase)
