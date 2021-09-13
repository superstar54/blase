# add
from blase.batoms import Batom
h1 = Batom('h2o', 'H_1', [[0, 0, 0], [2, 0, 0]])
h2 = Batom('h2o', 'H_2', [[0, 0, 2], [2, 0, 2]])
h = h1 + h2


from ase.build import molecule, fcc111
from blase.batoms import Batoms
h2o = Batoms(atoms=molecule('H2O'), label = 'h2o')
h2o.model_type = '1'


from blase.batoms import Batoms
from blase.batom import Batom
h2 = Batom('h2o', 'H', [[0, 0, 0], [1, 0, 0]])
h_new = h2.copy('h2o_new', 'H')


from blase.batoms import Batom
h = Batom('h2o', 'H', [[0, 0, 0], [2, 0, 0]])
h1 = h.copy('h2o', 'H_1')
h1.translate([0, 0, 2])
h.extend(h1)

from ase.build import molecule, fcc111
from blase.batoms import Batoms
pt111 = fcc111('Pt', (5, 5, 4), vacuum = 5.0)
pt111 = Batoms(atoms = pt111, label = 'pt111')
pt111.replace('Pt', 'Au', [93])
atoms = molecule('H2O')
h2o = Batoms(atoms = atoms, label = 'h2o')


# polyhedra
from blase.batoms import Batoms
from ase.io import read
atoms = read('docs/source/_static/datas/tio2.cif')
tio2 = Batoms(label = 'tio2', atoms = atoms, model_type = '2', polyhedra_dict = {'Ti': ['O']}, color_style="VESTA")



# materials type
from ase.build import molecule
from blase.batoms import Batoms
h2o = Batoms(label = 'h2o', atoms=molecule('H2O'), material_style='mirror')
h2o['H'][[0, 1]]
h2o['H'][[0, 1]] = [[3, 0, 0], [0, -3, 0]]

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
from blase.batom import Batom
c = Batom('co', 'C', [[0, 0, 0], [1.2, 0, 0]])
c.repeat([3, 3, 3], np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]))

# repeat
from ase.build import molecule, fcc111
from blase.batoms import Batoms
import numpy as np
pt111 = fcc111('Pt', (1, 1, 4), vacuum = 5.0)
pt111.pbc = True
pt111 = Batoms(atoms = pt111, label = 'pt111')
pt111 *= [4, 4, 1]
atoms = molecule('CO')
co = Batoms(atoms = atoms, label = 'co')
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
pt111 = Batoms(atoms = pt111, label = 'pt111')
pt111 *= [4, 4, 1]
pt111.set_cell(pt111.get_cell()*1.2, scale_atoms=True)


# replace
from ase.build import molecule, fcc111
from blase.batoms import Batoms
import numpy as np
pt111 = fcc111('Pt', (1, 1, 4), vacuum = 5.0)
pt111.pbc = True
pt111 = Batoms(atoms = pt111, label = 'pt111')
pt111.replace('Pt', 'Au', [2, 3])


@ index
from ase.build import molecule, fcc111
from blase.batoms import Batoms
import numpy as np
pt111 = fcc111('Pt', (1, 1, 4), vacuum = 5.0)
pt111.pbc = True
pt111 = Batoms(atoms = pt111, label = 'pt111')

print(pt111['Pt'][0])
pt111['Pt'][0] = [0, 0, 2.0]

index = [2, 3]
print(pt111['Pt'][index])
pt111['Pt'][index][:, 0] += 2


