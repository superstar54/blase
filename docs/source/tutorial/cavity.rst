
============================
Cavity
============================

The :meth:`~Batoms.draw_cavity` function is to draw cavity in porous materials. We find the cavity with a grid based algorithm. Here we search cavity with radius larger than 9.0 in the MOF-5 crystal.

>>> from ase.io import read
>>> from blase.batoms import Batoms
>>> atoms = read('docs/source/_static/datas/mof-5.cif')
>>> mof = Batoms(label = 'mof-5', atoms = atoms)
>>> mof.draw_cell()
>>> mof.model_type = 2
>>> mof.draw_cavity(9.0)
>>> mof.render.light_energy = 3
>>> mof.render.run([1, 0, 0], engine = 'eevee')


.. image:: ../_static/cavity.png
   :width: 8cm


