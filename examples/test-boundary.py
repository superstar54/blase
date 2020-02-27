from ase.build import bulk
from blase.tools import write_blender

atoms = bulk('Pt', cubic = True)
atoms = atoms*[6, 6, 6]

kwargs = {'show_unit_cell': 1, 
          'engine': 'BLENDER_WORKBENCH', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH', CYCLES
          'radii': 1.0,
          # 'display': True,
          'boundary_list': [{'d': 10.0, 'index': [1, 1, 0]}],
          'outfile': 'figs/test-boundary',
          }
write_blender(atoms, **kwargs)
