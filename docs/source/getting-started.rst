**********************
Getting Started
**********************

There are three ways to access the functionality of blase:

Option 1: Use blaseio to run Blender in the background
===========================================

>>> from ase.build import molecule
>>> from blaseio import write_blender
>>> atoms = molecule('H2O')
>>> batoms = {'atoms': atoms, 'model_type': '1'}
>>> blase = {'output_image': 'figs/h2o',}
>>> write_blender(batoms, blase)

|h2o|




Option 2: Import blase as Python module inside Blender’s Python console
===========================================


|h2o-2|


Option 3: Use the blase's GUI
===========================================

|h2o-3|


.. |h2o| image:: ../_static/h2o.png
   :width: 3cm
.. |h2o-2| image:: ../_static/h2o-2.png
   :width: 20cm
.. |h2o-3| image:: ../_static/h2o-3.png
   :width: 20cm