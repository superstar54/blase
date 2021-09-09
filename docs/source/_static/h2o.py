from ase.build import molecule
from blaseio import write_blender
atoms = molecule('H2O')
batoms = {'atoms': atoms, 'model_type': '1'}
blase = {'output_image': 'figs/h2o',}
write_blender(batoms, blase)

