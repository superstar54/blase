**********************
Getting Started
**********************


Start Blender, and open a Python console, and run

>>> from ase.build import molecule
>>> from blase.batoms import Batoms
>>> atoms = molecule('H2O')
>>> h2o = Batoms(label = 'h2o', atoms = atoms)


.. image:: _static/h2o-2.png
   :width: 20cm

Rendering the image

>>> h2o.render.run(direction = [1, 0 ,0], engine = 'eevee')

.. image:: _static/h2o.png
   :width: 5cm

Use blaseio to run Blender in the background
==============================================

>>> from ase.build import molecule
>>> from blaseio import write_blender
>>> atoms = molecule('H2O')
>>> batoms = {'atoms': atoms}
>>> blase = {'output': 'h2o.png',}
>>> write_blender(batoms, blase)


