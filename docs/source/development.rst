.. _devel:

============
Development
============

In Blender, when the undo/redo operation system (ctrl+z / ctrl+shift+z) is performed, all objects in the scene are fully recreated. Pointers to the old objects will fail. Therefore, please store list of object names instead of object itself, and then use them by ``bpy.data.objects[name]``. This is same for collection and others.


IDE
=======

``VS Code`` with ``Blender Development`` extenstion is highly recommended to used . Please read : https://marketplace.visualstudio.com/items?itemName=JacquesLucke.blender-development

Known bugs
===================




Development topics
=====================



Polyhedra
----------------
  
- Information about coordination polyhedra: volumes
  
Atom
-----------

- Add vectors (arrows) to atoms to represent magnetic moments and so on
- boundary cut the original atoms in batoms, how to handle?
- Add particle system

Light
----------

- For perspective view, ``SUN`` light bound to camera is not good.
  

Algorithm
------------------

- Find high order bond for aromatics. (eg. using CSD Python API, https://downloads.ccdc.cam.ac.uk/documentation/API/descriptive_docs/substructure_searching.html#)
- Add 2D slices of volumetric data in their 3D image.
- Faster bond-search algorithm
- Cavity