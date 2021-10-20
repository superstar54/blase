.. blase documentation master file, created by
   sphinx-quickstart on Wed Sep  8 10:33:04 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=================================
Welcome to blase's documentation!
=================================
Blase is a Python package for editing and rendering atoms and molecules objects using blender. A Python interface that allows for automating workflows.

>>> from blase.batoms import Batoms
>>> h2o = Batoms('h2o', {'O': [[0, 0, 0.40]], 'H': [[0, -0.76, -0.2], [0, 0.76, -0.2]]})

|logo|

Features:


* Model: space-filling, ball-stick, polyhedral, cavity and so on.
* File type: cif, xyz, cube, pdb, json, VASP-out and so on.
* Volumetric data (Isosurface)
* Animation
* GUI
* ``Flexible``: Python script, run interactively or in background.
* ``High quality rendering``:  3D models
* ``Free, Open Source``: Easy to download and install.
* ``Cross-platform``: (Linux, Windows, macOS)


.. toctree::
   :maxdepth: 3
   
   install
   getting-started
   tutorial/index
   module/index
   advanced/index
   gui/index
   gallery
   tips
   faqs
   development



* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _feedback: 
.. _affiliated packages: 

.. |logo|  image:: _static/batoms-h2o.png
   :width: 3cm