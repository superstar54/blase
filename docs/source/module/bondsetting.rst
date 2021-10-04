.. module:: blase.bondsetting

========================
The Bondsetting object
========================

By defaut, we use `ase.neighborlist.neighbor_list` to get the bond paris. We use default radius for every atoms, and if two spheres overlap, atoms are connected by a bond.

>>> from blase.bio import read
>>> tio2 = read('docs/source/_static/datas/tio2.cif')
>>> tio2.model_type = 1

.. image:: ../_static/bond-tio2.png
   :width: 8cm

You can print the default bondsetting by:

>>> tio2.bondsetting

.. image:: ../_static/bondsetting-tio2.png
   :width: 15cm

To build up coordination polyhedra, the value for ``polyhedra`` should be set to ``True``. To change setting for a bond pair by:

>>> tio2.bondsetting[('Ti', 'O')] = [0.5, 2.5, True, False]
>>> tio2.model_type = 2

.. image:: ../_static/polyhedra-tio2.png
   :width: 8cm


Search bond mode
==================
 
 * Do not search atoms beyond the boundary
  

>>> tio2.boundary = 0.01
>>> tio2.model_type = 2

.. image:: ../_static/search-bond-tio2-1.png
   :width: 8cm


* Search additional atoms if species1 is included in the boundary, the value for ``Search_bond`` should be set to ``True``. To change setting for a bond pair by:

>>> tio2.bondsetting[('Ti', 'O')] = [0.5, 2.5, True, True]
>>> tio2.update_boundary()
>>> tio2.model_type = 2


.. image:: ../_static/search-bond-tio2-2.png
   :width: 8cm


List of all Methods
===================

.. autoclass:: Bondsetting
   :members: