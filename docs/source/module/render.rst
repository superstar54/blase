.. module:: blase.render

===================
The Render object
===================

The :class:`Render` object is to render atomic structure.

>>> from ase.build import molecule
>>> from blase.batoms import Batoms
>>> atoms = molecule('C2H6SO')
>>> c2h6so = Batoms(label = 'c2h6so', atoms = atoms)
>>> c2h6so.render.run(direction = [1, 0, 0], engine = 'eevee')


Here, the first keyword ``direction`` specifies the direction of camera. Keyword ``engine`` specifies the render engine, and we used
the ``output_image`` keywords to specify name of output image.  Other
possible keywords are: ``canvas``, ``ratio``, ``camera_loc``, ``camera_type``, ``camera_lens``,
``ortho_scale``, ``studiolight``, ``light_loc``, ``light_type``, ``light_energy``, ``resolution_x`` and ``use_motion_blur``.

You can print the default bondsetting by:

>>> c2h6so.rendersetting()


Parameters
------------

- direction
- resolution_x
- ortho_scale
- ratio

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
   * -  .. image:: ../_static/render_workbench.png 
     -  .. image:: ../_static/render_eevee.png 
     -  .. image:: ../_static/render_cycles.png 


Viewpoint
===========
The ``direction`` keyword is used to set the direction of camera.


.. list-table::
   :widths: 25 25 25 25

   * - Top
     - Front
     - Right
     - Other direction
   * - [0, 0, 1]
     - [1, 0, 0]
     - [0, 1, 0]
     - [1, -0.3, 0.1]
   * -  .. image:: ../_static/render_direction_top.png 
     -  .. image:: ../_static/render_direction_front.png 
     -  .. image:: ../_static/render_direction_right.png 
     -  .. image:: ../_static/render_direction_any.png 



Edit the studio light
=======================

The Workbench engine use the studio light instead of the lights in the scene. You can choose the following studio light



.. list-table::
   :widths: 25 25 25 25 25 25

   * - ``Default``
     - ``basic.sl``
     - ``outdoor.sl``
     - ``paint.sl``
     - ``rim.sl``
     - ``studio.sl``
   * -  .. image:: ../_static/studiolight_Default.png 
     -  .. image:: ../_static/studiolight_basic.sl.png 
     -  .. image:: ../_static/studiolight_outdoor.sl.png 
     -  .. image:: ../_static/studiolight_paint.sl.png 
     -  .. image:: ../_static/studiolight_rim.sl.png 
     -  .. image:: ../_static/studiolight_studio.sl.png 


>>> h2o.render.run(studiolight = 'paint.sl')

To make your own studio light. Please read: https://docs.blender.org/manual/en/latest/editors/preferences/lights.html#prefs-lights-studio. Please save it, for example ``myown.sl``, then use it by:

>>> h2o.render.run(studiolight = 'mywon.sl')




Other methods
=============

* :meth:`~Render.set_camera`
* :meth:`~Render.set_light`
* :meth:`~Render.draw_plane`
* :meth:`~Render.export`
* :meth:`~Render.render`



List of all Methods
===================

.. autoclass:: Render
   :members: