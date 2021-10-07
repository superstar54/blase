.. _render:

====================
Rendering structure
====================

After you have created a structure, you can export nice images for publications or presentations by rendering.

The :mod:`Render <blase.render>` object controls various settings such as the viewpoint and the resolution of the generated image.


>>> from ase.build import molecule
>>> from blase.batoms import Batoms
>>> atoms = molecule('C2H6SO')
>>> c2h6so = Batoms(label = 'c2h6so', atoms = atoms)
>>> c2h6so.render.run(direction = [1, 0, 0], engine = 'eevee', resolution_x = 1000, output = 'c2h6so.png')

Possible keywords:

- ``output``, name of image.
- ``resolution_x``, resolution of image, default 1000.
- ``light_energy`` keyword to make the image lighter or darker, default 5.0.


Viewpoint
=============
The ``direction`` keyword is used to set the direction of camera.


.. list-table::
   :widths: 25 25 25 25

   * - Top view (``z``)
     - Front view (``x``)
     - Right view (``y``)
     - Other direction (``x, y, z``)
   * - [0, 0, 1]
     - [1, 0, 0]
     - [0, 1, 0]
     - [1, -0.3, 0.1]
   * -  .. image:: ../../_static/render_direction_top.png 
     -  .. image:: ../../_static/render_direction_front.png 
     -  .. image:: ../../_static/render_direction_right.png 
     -  .. image:: ../../_static/render_direction_any.png 


Engine
=============

Three engines can be used:

.. list-table::
   :widths: 25 25 25

   * - BLENDER_WORKBENCH
     - BLENDER_EEVEE
     - CYCLES
   * - Fast
     - Standard
     - Most powerfull, slow
   * -  .. image:: ../../_static/render_workbench.png 
     -  .. image:: ../../_static/render_eevee.png 
     -  .. image:: ../../_static/render_cycles.png 

.. toctree::
   :maxdepth: 1
   
   camera
   light
