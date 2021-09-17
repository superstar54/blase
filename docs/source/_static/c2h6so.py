from ase.build import molecule
from blaseio import write_blender

atoms = molecule('C2H6SO')
atoms.write('datas/c2h6so.xyz')
batoms = {'label': 'c2h6so', 'atoms': atoms, 'model_type': '1'}
blase = {'output_image': 'figs/c2h6so',}
write_blender(batoms, blase, display=True)
