================================
Volumetric data and isosurface
================================

The :mod:`Isosurfacesetting <blase.isosurfacesetting>` object is used to store and set all parameters related with volumetric data. It should always bind with a :class:`Batoms` object. Possible keywords are: ``level``, ``color``. 


Here we show a example of draw isosurfaces from cube file.

>>> from ase.io.cube import read_cube_data
>>> from blase.batoms import Batoms
>>> volume, atoms = read_cube_data('docs/source/_static/datas/h2o-homo.cube')
>>> h2o = Batoms('h2o', atoms = atoms, volume = volume, draw = False)


You can print the default isosurfacesetting by:

>>> h2o.isosurfacesetting

One add a isosurfacesetting by:

>>> h2o.isosurfacesetting[2] = [-0.002, [0, 0, 0.8, 0.5]]
>>> h2o.draw_isosurface()

The first value ``-0.002`` is the level for the isosurface, the second value ``[0, 0, 0.8, 0.5]`` is the color. The last value of the color ``0.5`` is used to set the ``transparency``.

.. image:: ../_static/volume_h2o.png
   :width: 10cm

