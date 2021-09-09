.. module:: blase.batoms

===================
The Batoms object
===================

The :class:`Batoms` object is a collection of objects. Here is how to define a H2O molecule.

>>> from ase.build import molecule
>>> from blase.batoms import Batoms
>>> atoms = molecule('H2O')
>>> h2o = Batoms(atoms = atoms, name = 'h2o', model_type = '1')
>>> h2o.draw()
>>> h2o.render()

.. image:: _static/h2o.png
   :width: 3cm

Here, the first keyword ``atoms`` specifies the ase ``Atoms`` object, and we used
the ``model_type`` keywords to specify model type.  Other
possible keywords are: ``boundary``, ``show_unit_cell``, ``isosurface``, ``kind_props``,
``color`` and ``draw``.


model_type
===================

Here, four models can be set:

- 0: Space-filling

- 1: Ball-and-stick

- 2: Polyhedral

- 3: Stick

>>> h2o = Batoms(atoms = atoms, model_type = '1')

or,

>>> h2o.set_model_type('1')


Other methods
=============

* :meth:`~Batoms.translate`
For example, move h2o molecule by a vector [0, 0, 5],

>>> h2o = translate([0, 0, 5])

* :meth:`~Batoms.rotate`

For example, rotate h2o molecule 90 degree around 'Z' axis:

>>> h2o.rotate(90, 'Z')

* :meth:`~Batoms.copy`
  
For example, copy h2o molecule:
        
>>> h2o_new = h2o.copy(name = 'h2o_new')

* :meth:`~Batoms.delete`

For example, delete the second atom in h2o molecule. Please note that index start from 0.

>>> h2o.delete([1])

* :meth:`~Batoms.replace`

For example, replace the all H in h2o molecule by S.

>>> h2o.replace([1, 2], 'S')

* :meth:`~Batoms.repeat`

>>> from ase.build import bulk
>>> from blase.batoms import Batoms
>>> au = bulk('Au', cubic = True)
>>> au = Batoms(atoms = au)
>>> au.draw()
>>> au.repeat([2, 2, 2])


* :meth:`~Batoms.extend`

>>> from ase.build import molecule, fcc111
>>> from blase.batoms import Batoms
>>> import numpy as np
>>> co = molecule('CO')
>>> co = Batoms(atoms = co, draw = True)
>>> au = fcc111('Au', (5, 5, 4), vacuum=5.0)
>>> au = Batoms(atoms = au, draw = True)
>>> co.translate(au.atoms[-1].position + np.array([0, 0, 2]))
>>> au.extend(co)

or,

>>> au = au + co


* :meth:`~Batoms.write`

Save atoms to file, please vist write method in ASE, https://wiki.fysik.dtu.dk/ase/ase/io/io.html?highlight=write#ase.io.write

>>> au.write('au111-co.cif')


* :meth:`~Batoms.render`

Render the atoms, and save to a png image.

>>> h2o.render(resolution_x = 1000, output_image = 'h2o.png')


List of all Methods
===================

.. autoclass:: Batoms
   :members: