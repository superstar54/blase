from ase.io import read
from blase.batoms import Batoms
atoms = read('docs/source/_static/datas/mof-5.cif')
mof = Batoms(label = 'mof-5', atoms = atoms)
mof.draw_cavity(9.0)
mof.draw_cell()
mof.model_type = 2
mof.render.light_energy = 3
mof.render.run([1, 0, 0], engine = 'eevee')
