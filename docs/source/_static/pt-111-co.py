from ase.build import fcc111, add_adsorbate
from ase.atoms import Atoms
from ase.build import bulk
from runblase import write_blender

bulk = bulk('Pt')
bulk.write('pt.in')
adsorbate = Atoms('CO')
adsorbate[1].z = 1.1
atoms = fcc111('Pt', (2, 2, 3), a=3.96, vacuum=7.0)
add_adsorbate(atoms, adsorbate, 1.8, 'ontop')
# atoms = atoms*[5, 5, 1]

batoms = {
        'model_type': '1',
        'kind_props' : {'Pt':{'scale': 1.0}, 'C':{'scale': 0.6}, 'O': {'scale': 0.6}},
        'show_unit_cell': True,
        'remove_bonds': {'Pt':['Pt']},
        }
blase = {
          'engine': 'BLENDER_WORKBENCH', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH', CYCLES
          'output_image': 'figs/pt-111-co',
  }

write_blender(atoms, batoms, blase, display = True)
