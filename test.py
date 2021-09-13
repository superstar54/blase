from numpy.core.defchararray import index
from ase.build import molecule, fcc111
from blase.batoms import Batoms
h2o = Batoms(atoms=molecule('H2O'), name = 'h2o')


from blase.batoms import Batoms
from blase.batom import Batom
h2 = Batom('h2o', 'H', [[0, 0, 0], [1, 0, 0]])
h_new = h2.copy('h2o_new', 'H')


from blase.batoms import Batom
h = Batom('h2o', 'H', [[0, 0, 0], [2, 0, 0]])
h1 = h.copy('h2o', 'H_1')
h1.translate([0, 0, 2])
h.extend(h1)

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
co = Batoms(atoms = co)
co.translate([0, 0, 2])
h2o = molecule('H2O')
h2o = Batoms(atoms = h2o)
h2o.extend(co)


# repeat
from ase.build import molecule, fcc111
from blase.batoms import Batoms
import numpy as np
pt111 = fcc111('Pt', (1, 1, 4), vacuum = 5.0)
pt111.pbc = True
pt111 = Batoms(atoms = pt111, name = 'pt111')
pt111 *= [4, 4, 1]
atoms = molecule('CO')
co = Batoms(atoms = atoms, name = 'co')
t = pt111['Pt'].positions[0] - co['C'].positions[0] + np.array([0, 0, 2])
co.translate(t)
pt111.extend(co)
pt111.write('pt111-co.in')

# set_cell
from ase.build import molecule, fcc111
from blase.batoms import Batoms
import numpy as np
pt111 = fcc111('Pt', (1, 1, 4), vacuum = 5.0)
pt111.pbc = True
pt111 = Batoms(atoms = pt111, name = 'pt111')
pt111 *= [4, 4, 1]
pt111.set_cell(pt111.get_cell()*1.2, scale_atoms=True)


# replace
from ase.build import molecule, fcc111
from blase.batoms import Batoms
import numpy as np
pt111 = fcc111('Pt', (1, 1, 4), vacuum = 5.0)
pt111.pbc = True
pt111 = Batoms(atoms = pt111, name = 'pt111')
pt111.replace('Pt', 'Au', [2, 3])


@ index
from ase.build import molecule, fcc111
from blase.batoms import Batoms
import numpy as np
pt111 = fcc111('Pt', (1, 1, 4), vacuum = 5.0)
pt111.pbc = True
pt111 = Batoms(atoms = pt111, name = 'pt111')

print(pt111['Pt'][0])
pt111['Pt'][0] = [0, 0, 2.0]

index = [2, 3]
print(pt111['Pt'][index])
pt111['Pt'][index][:, 0] += 2


