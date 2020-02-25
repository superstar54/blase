from ase.cluster import wulff_construction
from ase.visualize import view
from blase.tools import get_bondpairs, write_blender

surfaces = [(1, 1, 1), (1, 0, 0)]
energies = [1.28, 1.69]
atoms = wulff_construction('Au', surfaces, energies, 5000, 'fcc')
# view(atoms1)
for i in range(3):
  atoms.positions[:, i] -= min(atoms.positions[:, i])
atoms.positions[:, i] -= min(atoms.positions[:, i])
atoms.center(vacuum=2.0)
# view(atoms)
camera_loc =  atoms.get_center_of_mass() + [0, -300, 0]


kwargs = {
          'engine': 'BLENDER_WORKBENCH', #'CYCLES', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH'
          'camera_loc': camera_loc,  # distance from camera to front atom
          'camera_type': 'PERSP',  #  ['PERSP', 'ORTHO', 'PANO']
          'camera_lens': 100,  #
          'camera_target': atoms.get_center_of_mass(), #
          'radii': 1.0, 
          'display': True,
          'outfile': 'figs/wulff'}   
write_blender(atoms, **kwargs)