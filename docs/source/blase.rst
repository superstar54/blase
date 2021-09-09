.. module:: blase.bio

===================
The Blase object
===================

The :class:`Blase` object is to render atomic structure.

>>> from ase.build import molecule
>>> from blase.batoms import Batoms
>>> from blase.bio import Blase
>>> atoms = molecule('H2O')
>>> h2o = Batoms(atoms = atoms, model_type = '1')
>>> h2o.draw()
>>> bobj = Blase(engine = 'BLENDER_WORKBENCH',
                 output_image = 'h2o.png',
                 bbox = 'all')
>>> bobj.render()


Here, the first keyword ``engine`` specifies the render engine, and we used
the ``output_image`` keywords to specify name of output image.  Other
possible keywords are: ``bbox``, ``camera``, ``camera_loc``, ``camera_type``, ``camera_lens``,
``light_loc``, ``light_type`` and ``resolution_x``.


engine
=============

Here, three models can be set:

- BLENDER_WORKBENCH
- BLENDER_EEVEE
- CYCLES


Other methods
=============

* :meth:`~Blase.set_camera`
* :meth:`~Blase.set_light`
* :meth:`~Blase.draw_plane`
* :meth:`~Blase.export`
* :meth:`~Blase.render`



List of all Methods
===================

.. autoclass:: Blase
   :members: