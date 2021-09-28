.. module:: blase.batoms

===================
How to load Batoms from a collection?
===================

When we restart from a ``.blend`` file. We can rebuild the :class:`Batoms` objects from collections. Keyword ``from_collection`` in :class:`Batoms` object is uesd to load a :class:`Batoms` collection by the name.


Molecule:

>>> from blase.batoms import Batoms
>>> co = Batoms(from_collection = 'co')

.. image:: _static/ase-co.png
   :width: 3cm

Crystal:

>>> from blase.batoms import Batoms
>>> fe = Batoms(from_collection = 'fe')

.. image:: _static/ase-fe.png
   :width: 3cm

.. module:: blase.batoms
