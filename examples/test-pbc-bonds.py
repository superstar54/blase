#!/usr/bin/env python
from ase.build import surface, fcc111
from ase.visualize import view
from blase.tools import get_bondpairs, write_blender
#============================================================
atoms = fcc111('Pt', (2, 2, 1), vacuum=5.0)
# atoms.pbc = [True, True, True]
# view(atoms)
camera_loc = atoms.get_center_of_mass() + [0, -60, 30]

kwargs = {'show_unit_cell': True, 
          'camera_loc': camera_loc,  # distance from camera to front atom
          'camera_type': 'PERSP',  #  ['PERSP', 'ORTHO', 'PANO']
          'camera_lens': 100,  #
          'camera_target': atoms.get_center_of_mass(), #
          'ortho_scale': None, #
          'radii': 0.6,
          'bond_cutoff': 1.0,
          'resolution_x': 1000,
          'display': True,
          'outfile': 'figs/test-pbc-bonds',
          }
write_blender(atoms, **kwargs)