
===================
Hydrogen bond
===================

By defaut, we use `ase.neighborlist.neighbor_list` to get the bond paris. We use default radius for every atoms, and if two spheres overlap, atoms are connected by a bond.

You can print the default bondsetting by:



To build up coordination polyhedra, the value for ``polyhedra`` should be set to ``True``. To change setting for a bond pair by:

>>> tio2.bondsetting[('Ti', 'O')] = [2.5, True, False]


 Search mode:
 
 * Search species2 bonded to species1



