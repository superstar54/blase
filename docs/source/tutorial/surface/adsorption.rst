
==========================================
Surface adsorption
==========================================

In this tutorial we will adsorb a CO molecule on Au(111) surface.

>>> import numpy as np
>>> from ase.build import fcc111, molecule
>>> from blase.batoms import Batoms
>>> atoms = fcc111('Au', size = (4, 4, 4), vacuum=0)
>>> au111 = Batoms(label = 'au111', atoms = atoms)
>>> au111.cell[2, 2] += 10
>>> co = Batoms(label = 'co', atoms = molecule('CO'))
>>> # In Edit mode, you can see the index of gold atoms.
>>> t = au111['Au'][58] - co['C'][0] + np.array([0, 0, 2])
>>> co.translate(t)
>>> au111 = au111 + co

Save structure to file, POSCAR, xyz and so on.

>>> au111.write('POSCAR')

Save image, as show bellow:

>>> au111.render([1, 0.1, 0.1], engine = 'eevee')


.. image:: ../../_static/au111-co.png
   :width: 8cm


