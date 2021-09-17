.. _devel:

===========
Development
===========

In Blender, when the undo/redo operation system (ctrl+z / ctrl+shift+z) is performed, all objects in the scene are fully recreated. Pointers to the old objects will fail. Therefore, please store list of object names instead of object itself, and then use them by ``bpy.data.objects[name]``. This is same for collection and others.



Known bugs:




Development topics:



- Add vectors (arrows) to atoms to represent magnetic moments and so on
- measure: interatomic distance, angle and dihedral angles
- information about coordination polyhedra: volumes



Algorithm:

- add 2D slices of volumetric data in their 3D image.
- Bond-search algorithm

