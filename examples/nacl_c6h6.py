import numpy as np
from ase import Atoms
from ase.visualize import view
from ase.build import molecule
from blase.tools import get_bondpairs, write_blender

a = 5.64  # Lattice constant for NaCl
cell = [a / np.sqrt(2), a / np.sqrt(2), a]
atoms = Atoms(symbols='Na2Cl2', pbc=True, cell=cell,
              scaled_positions=[(.0, .0, .0),
                                (.5, .5, .5),
                                (.5, .5, .0),
                                (.0, .0, .5)]) * (3, 4, 2) + molecule('C6H6')

# Move molecule to 3.5Ang from surface, and translate one unit cell in xy
atoms.positions[-12:, 2] += atoms.positions[:-12, 2].max() + 3.5
atoms.positions[-12:, :2] += cell[:2]
atoms.cell[0][0] += 10
# Mark a single unit cell
# atoms.cell = cell
# view(atoms)

# View used to start ag, and find desired viewing angle
#view(atoms)
# rot = '35x,63y,36z'  # found using ag: 'view -> rotate'


atoms.rotate(90, 'y')
# view(atoms)
bondatoms = get_bondpairs(atoms, cutoff=1.2, rmbonds=[['Na', 'Na'], ['Cl', 'Cl']])
# print(bondatoms)
kwargs = {
          # 'show_unit_cell': 1, 
          'radii': 0.6, 
          # 'functions': [['draw_plane', {'size': 200, 'loc': (0, 0, -1.0)}]],
          'bonds': 'all',
          'display': True,    
          'radii': {'C': 0.8, 'H': 0.8}, 
          'bondlist': bondatoms,
          'bondlinewidth': 0.4,  # radius of the cylinders representing bonds
          'outfile': 'figs/NaCl_C6H6',
          }
write_blender(atoms, **kwargs)
