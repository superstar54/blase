import numpy as np
from ase import Atoms
from ase.visualize import view
from ase.build import molecule
from runblase import write_blender

a = 5.64  # Lattice constant for NaCl
cell = [a / np.sqrt(2), a / np.sqrt(2), a]
nacl = Atoms(symbols='Na2Cl2', pbc=True, cell=cell,
              scaled_positions=[(.0, .0, .0),
                                (.5, .5, .5),
                                (.5, .5, .0),
                                (.0, .0, .5)]) * (3, 4, 2)
c6h6 = molecule('C6H6')

# Move molecule to 3.5Ang from surface, and translate one unit cell in xy
c6h6.positions[:, 1] += nacl.positions[:, 1].max() + 3.5
c6h6.positions[:, :1] += cell[:1]
nacl.cell[1][1] += 0


batoms1 = {'atoms': nacl,
        'model_type': '0',
        'show_unit_cell': True,
        }
batoms2 = {'atoms': c6h6,
        'model_type': '1',
        'show_unit_cell': True,
        }
blase = {
          'engine': 'BLENDER_WORKBENCH', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH', CYCLES
          'output_image': 'figs/NaCl_C6H6',
  }

write_blender([batoms1, batoms2], blase)
