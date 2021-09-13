.. module:: blase.batom

===================
The Batom object
===================

The :class:`Batom` object is a object for one species. Here is how to define a H species

>>> from blase.batom import Batom
>>> h = Batom(label = 'h2o', species = 'H', positions = [[0, 0, 0], [2.0, 0. 0]])

.. image:: _static/h.png
   :width: 3cm

Here, the ``label`` keywords to specify the name, and ``species`` keywords 
to specify the species ``H``, and the ``positions`` keywords to 
specify their positions of the H atoms. Other possible keywords are: ``element``, ``scale``, ``kind_props``, ``color_style``,
``material_style``, ``bsdf_inputs`` and ``draw``.




Other methods
=============

* :meth:`~Batom.translate`
For example, move all ``H`` species by a vector [0, 0, 5],

>>> h = translate([0, 0, 5])

* :meth:`~Batom.rotate`

For example, rotate ``H`` species 90 degree around 'Z' axis:

>>> h.rotate(90, 'Z')

* :meth:`~Batom.copy`
  
For example, copy ``H`` species:
        
>>> h_new = h.copy('h_new', 'H')

* :meth:`~Batom.delete`

For example, delete the second and the third atom in ``H`` species. Please note that index start from 0.

>>> h.delete([1, 2])

Or,

>>> del h[[1, 2]]

* :meth:`~Batom.replace`

For example, replace the all H in h molecule by S.

>>> h.replace('H', 'S', [1])

* :meth:`~Batom.repeat`

>>> from ase.build import bulk
>>> from blase.batom import Batom
>>> au = bulk('Au', cubic = True)
>>> au = Batom(atoms = au)
>>> au.draw()
>>> au.repeat([2, 2, 2])


* :meth:`~Batom.extend`

>>> from ase.build import molecule, fcc111
>>> from blase.batom import Batom
>>> import numpy as np
>>> co = molecule('CO')
>>> co = Batom(atoms = co, draw = True)
>>> au = fcc111('Au', (5, 5, 4), vacuum=5.0)
>>> au = Batom(atoms = au, draw = True)
>>> co.translate(au.atoms[-1].position + np.array([0, 0, 2]))
>>> au.extend(co)

or,

>>> au = au + co


* :meth:`~Batom.write`

Save atoms to file, please vist write method in ASE, https://wiki.fysik.dtu.dk/ase/ase/io/io.html?highlight=write#ase.io.write

>>> au.write('au111-co.cif')

* :meth:`~Batom.show_index`

Show the index of atoms.

>>> au.show_index(index_type = 0)

* :meth:`~Batom.render`

Render the atoms, and save to a png image.

>>> h2o.render(resolution_x = 1000, output_image = 'h2o.png')


List of all Methods
===================

.. autoclass:: Batom
   :members: