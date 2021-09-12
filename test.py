from ase.build import molecule, fcc111
from blase.batoms import Batoms
atoms = molecule('H2O')
h2o = Batoms(atoms = atoms, name = 'h2o')
h2o2 = h2o.copy('h2o2')


pt111 = fcc111('Pt', (5, 5, 4), vacuum = 5.0)
pt111 = Batoms(atoms = pt111, name = 'pt111')
pt111.replace('Pt', 'Au', [93])
atoms = molecule('H2O')
h2o = Batoms(atoms = atoms, name = 'h2o')


# from blase.batoms import Batoms
# from ase.io import read
# atoms = read('docs/source/_static/datas/tio2.cif')
# tio2 = Batoms(atoms = atoms, name = 'tio2', model_type = '2', polyhedra_dict = {'Ti': ['O']}, color="VESTA")



# extend
from ase.build import molecule, fcc111
from blase.batoms import Batoms
import numpy as np
co = molecule('CO')
co = Batoms(atoms = co, draw = True)
co.translate([0, 0, 2])
h2o = molecule('H2O')
h2o = Batoms(atoms = h2o, draw = True)
h2o.extend(co)


# repeat
from ase.build import molecule, fcc111
from blase.batoms import Batoms
import numpy as np
pt111 = fcc111('Pt', (1, 1, 4), vacuum = 5.0)
pt111.pbc = True
pt111 = Batoms(atoms = pt111, name = 'pt111')
pt111.repeat([5, 5, 1])
atoms = molecule('CO')
co = Batoms(atoms = atoms, name = 'co')
t = pt111['Pt'].positions[31] - co['C'].positions[0] + np.array([0, 0, 2])
co.translate(t)
pt111.extend(co)
pt111.write('pt111-co.in')

