.. module:: blase.isosurfacesetting

=============================
The Isosurfacesetting object
=============================

The :class:`Isosurfacesetting` object is used to store and set all parameters related with volumetric data. It should always bind with a :class:`Batoms` object. Possible keywords are: ``level``, ``color``. 


Here we show a example of draw isosurfaces from cube file.

>>> from ase.io.cube import read_cube_data
>>> from blase.batoms import Batoms
>>> volume, atoms = read_cube_data('docs/source/_static/datas/h2o-homo.cube')
>>> h2o = Batoms('h2o', atoms = atoms, volume = volume, draw = False)


You can print the default isosurfacesetting by:

>>> h2o.isosurfacesetting

One add a isosurfacesetting by. 

>>> h2o.isosurfacesetting[1] = [0.002, [1, 1, 0, 0.5]]
>>> h2o.isosurfacesetting[2] = [-0.002, [0, 0, 0.8, 0.5]]
>>> h2o.draw_isosurface()

.. image:: ../_static/volume_h2o.png
   :width: 5cm


List of all Methods
====================

.. autoclass:: Isosurfacesetting
   :members: