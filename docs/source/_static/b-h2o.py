from ase.build import molecule
from blase.batoms import Batoms
h2o = molecule('H2O')
h2o = Batoms(atoms = h2o, name = 'h2o', model_type = '1')
h2o.draw()
#
co = molecule('CO')
co = Batoms(atoms = co, name = 'co', model_type = '1')
co.draw()
#
h2o.translate([5, 0, 0])
co.rotate(90, 'X')

