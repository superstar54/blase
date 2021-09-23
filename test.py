#
from blase import Batoms
h2o = Batoms({'O': [[0, 0, 0.40]], 'H': [[0, -0.76, -0.2], [0, 0.76, -0.2]]})
co = Batoms({'C':[[0, 0, 0]], 'O':[[1.2, 0, 0]]})

from blase.cell import Bcell
cell = Bcell(label = 'pt', array = [2, 2, 2])

from blase import Batoms
a = 3.96
positions = [[0, 0, 0], [a/2, a/2, 0], [a/2, 0, a/2], [0, a/2, a/2]]
pt = Batoms({'Pt': positions}, pbc = True, cell = (a, a, a))
pt.translate([0, 0 ,5])
# pt.repeat([2, 2, 2])
pt1 = pt.copy('pt1')

# boundary
from blase import Batoms
a = 3.96
positions = [[0, 0, 0], [a/2, a/2, 0], [a/2, 0, a/2], [0, a/2, a/2]]
pt = Batoms({'Pt': positions}, pbc = True, cell = (a, a, a))
pt.boundary = 0.01


# polyhedra
from blase.batoms import Batoms
from blase.bio import read
tio2 = read('docs/source/_static/datas/tio2.cif')
tio2.boundary = 0.0
tio2.model_type = 1




tio2.bondsetting
tio2.model_type = 1
tio2.bondsetting[('Ti', 'O')] = [0.0, 2.5, True, False]
tio2.model_type = 2


# batom load_frames
from blase import Batom
import numpy as np
positions = np.array([[0, 0 ,0], [1.52, 0, 0]])
h = Batom('h2o', 'H', positions)
images = []
for i in range(10):
    images.append(positions + [i, 0, 0])

h.load_frames(images)

# batoms load_frames
from ase.io import read, write
from blase import Batoms
atoms = molecule('C2H6SO')
images = []
for i in range(20):
    temp = atoms.copy()
    temp.rotate(18*i, 'z')
    images.append(temp)

write('c2h6so-animation.xyz', images)
images = read('c2h6so-animation.xyz', index = ':')
c2h6so = Batoms(label = 'c2h6so', atoms = images)
c2h6so.load_frames()
c2h6so.render(animation = True)

# add
import batom
from blase.batoms import Batom
h1 = Batom('h2o', 'H_1', [[0, 0, 0], [1.52, 0, 0]])
h2 = Batom('h2o', 'H_2', [[0, 0, 2], [2, 0, 2]])
h = h1 + h2

h3=Batom(from_batom_object = 'atom_co_C')


# build from collection
from ase.build import molecule, fcc111
from blase import Batoms
atoms = molecule('H2O')
atoms.center(3.0)
atoms.pbc = True
atoms = atoms*[20, 20, 20]
h2o = Batoms(label = 'h2o', atoms = atoms)
h2o2 = Batoms(from_collection='h2o')
h = Batom(from_batom='atom_h2o_H')

# search boundary
from ase.build import bulk
from blase.batoms import Batoms
pt = bulk('Pt', cubic = True)
pt = Batoms(label = 'pt', atoms = pt)
pt.boundary =[-0.02, 1.02, -0.02, 1.02, -0.02, 1.02]

# delete
from ase.build import molecule
from blase import Batoms
h2o = Batoms(label = 'h2o', atoms = molecule('H2O'))
h2o_new = h2o.copy('h2o_new')
del h2o_new['H'][[0]]

#
# bondsetting
from ase.build import molecule
from blase import Batoms
c2h6so = Batoms(label = 'c2h6so', atoms = molecule('C2H6SO'))


from blase.bio import read
batoms = read('docs/source/_static/datas/perovskite.cif')
batoms.repeat([2, 2, 2])
batoms.bondsetting
batoms.bondsetting[('Pb', 'I')] = [0, 3.5, True, False]
batoms.boundary = 0.2
batoms.model_type = 2
batoms.render()

#cavity
from blase.bio import read
mof = read('docs/source/_static/datas/mof-5.cif')
mof.bondsetting[('Zn', 'O')] = [0, 2.5, True, False]
mof.bondsetting[('C', 'H')] = [0, 1.4, False, False]
mof.boundary = 0.01
mof.draw_cavity(9.0)
mof.model_type = 2
mof.render()

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





# get_distances
from blase.batoms import Batoms
from ase.io import read
atoms = read('docs/source/_static/datas/tio2.cif')
print(atoms.get_distances(0, [2 ,3]))
tio2 = Batoms(label = 'tio2', atoms = atoms, model_type = '2', polyhedra_dict = {'Ti': ['O']}, color_style="VESTA")
tio2.get_distances('Ti', 0, 'O', [0, 1])

# get_angles
from ase.build import molecule, fcc111
from blase.batoms import Batoms
atoms = molecule('H2O')
h2o = Batoms(atoms = atoms, label = 'h2o')
h2o.get_angle('H', 0, 'O', 0, 'H', 1)

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
co = Batoms(label = 'co', atoms = co)
co.translate([0, 0, 2])
h2o = molecule('H2O')
h2o = Batoms(label = 'h2o', atoms = h2o)
h2o.extend(co)


# repeat
from blase.batom import Batom
import numpy as np
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
pt111.set_cell(pt111.cell.array*1.2, scale_atoms=True)


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


from pymatgen.core.structure import Molecule
c_monox = Molecule(["C","O"], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.2]])
print(c_monox)

from blase import Batoms
from pymatgen.core import Lattice, Structure
fe = Structure(Lattice.cubic(2.8), ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])
fe = Batoms(label = 'fe', structure=fe)


# C60
from ase.build import molecule
from blase import Batoms
atoms = molecule('C60')
c60 = Batoms(label = 'C60', atoms = atoms)


from blase import Batoms
from ase.cluster.cubic import FaceCenteredCubic
surfaces = [(1, 0, 0), (1, 1, 1), (1, -1, 1)]
layers = [6, 5, -1]
atoms = FaceCenteredCubic('Cu', surfaces, layers)
nano = Batoms(label = 'nano', atoms = atoms)

import numpy as np

a = np.array([[0, 0.25, 0.5], [0, 0.5, 0.75], [0.25, 0.75, 0.5]])
np.where(a[:, 1]>0.3) and np.where(a[:, 2] > 0.5)
c = np.array([a-2, a-1, a, a + 1, a + 2, a + 3])

b = [-1.3, 3.3]

np.repeat()


# CO adsorption
import numpy as np
from ase.build import fcc111, molecule
from blase import Batoms
atoms = fcc111('Pt', size = (1, 1, 4), vacuum=0)
pt111 = Batoms(label = 'pt111', atoms = atoms)
pt111.repeat([4, 4, 1])
pt111.cell[2, 2] += 10
co = Batoms(label = 'co', atoms = molecule('CO'))
t = pt111['Pt'][27] - co['C'][0] + np.array([0, 0, 2])
co.translate(t)
pt111 = pt111 + co
pt111.write('pt111-co.in')
