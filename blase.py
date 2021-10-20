from blase.batoms import Batoms
from ase.build import fcc111
from ase.io import read
from blase.butils import removeAll
removeAll()

dbba = read('/home/xing/cp2k/surface/ag/nanographene/model/dbba--2h.xyz')
ag111 = fcc111('Ag', (1, 1, 4), vacuum = 0)

ag111 = Batoms('ag111', atoms = ag111)
ag111.cell[2, 2] += 10
dbba = Batoms('dbba', atoms = dbba)
dbba['C'].color = [1, 0, 0, 1]

ag = Batoms('ag', {'Ag_1': [[0, 0, 10]]})
ag['Ag_1'].color = [0, 1, 0, 1]

ag111.repeat([4, 5, 1])
ag111 = ag111 + ag + dbba
ag111.boundary = [2, 0, 0]
ag111.write('model/ag111-ag-dbba.in')