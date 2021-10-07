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

- Search atoms outside boundary




Development topics
=====================


Bond
------------
- add Hydrogen bond
- add different style for bond, dot and dash line, one color for all bonds.
  
Polyhedra
----------------
  
- add different style for Polyhedra, show polyhedra edges or not
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

- Add 2D slices of volumetric data in their 3D image.
- Faster bond-search algorithm
- Search additional atoms recursively, molecule crystal

