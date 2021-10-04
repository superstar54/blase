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
possible keywords are: ``canvas``, ``camera_loc``, ``camera_type``, ``camera_lens``,
``light_loc``, ``light_type``, ``light_energy`` and ``resolution_x``.

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



Direction
===========




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