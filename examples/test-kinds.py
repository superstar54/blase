#!/usr/bin/env python
from ase.build import surface, fcc111
from ase.visualize import view
from blase.tools import get_bondpairs, write_blender
#============================================================
atoms = fcc111('Pt', (7, 7, 3), vacuum=0.0)
kinds = []
kind_props = {
'Pt_0': {'color': [208/255.0, 208/255.0, 224/255.0]},
'Pt_1': {'color': [225/255.0, 128/255.0, 0/255.0]},
'Pt_2': {'color': [0/255.0, 191/255.0, 56/255.0]},
}
for i in range(len(atoms)):
     ind = int((atoms[i].x/5))
     kind = atoms[i].symbol + '_{0}'.format(ind)
     kinds.append(kind)
atoms.kinds = kinds

camera_loc = atoms.get_center_of_mass() + [0, -60, 30]

kwargs = {'show_unit_cell': 0, 
          'camera_loc': camera_loc,  # distance from camera to front atom
          'camera_type': 'PERSP',  #  ['PERSP', 'ORTHO', 'PANO']
          'camera_lens': 100,  #
          'camera_target': atoms.get_center_of_mass(), #
          'ortho_scale': None, #
          'kind_props': kind_props,
          'outfile': 'figs/test-kinds',
          }
write_blender(atoms, **kwargs)