
============================
Build structure for surface adsorption
============================

In this tutorial we will adsorb a CO molecule on Pt (111) surfaces.

>>> import numpy as np
>>> from ase.build import fcc111, molecule
>>> from blase import Batoms
>>> 
>>> atoms = fcc111('Pt', size = (1, 1, 4), vacuum=0)
>>> pt111 = Batoms(label = 'pt111', atoms = atoms)
>>> pt111.repeat([4, 4, 1])
>>> pt111.cell[2, 2] += 10
>>> co = Batoms(label = 'co', atoms = molecule('CO'))
>>> # In Edit mode, you can see the index of platinum atoms.
>>> t = pt111['Pt'][27] - co['C'][0] + np.array([0, 0, 2])
>>> co.translate(t)
>>> pt111 = pt111 + co
>>> # save atoms to file, POSCAR, xyz and so on.
>>> pt111.write('pt111-co.in')


.. image:: _static/pt111-co.png
   :width: 8cm


