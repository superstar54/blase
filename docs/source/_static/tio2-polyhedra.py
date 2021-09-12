from ase.io import read
from blaseio import write_blender

atoms = read('datas/tio2.cif')
batoms = {'atoms': atoms,
        'model_type': '2',
        'polyhedra_dict': {'Ti': ['O']},
        'color': 'VESTA',
        }
blase = {
          'output_image': 'tio2-polyhedra',
  }
write_blender(batoms, blase)
