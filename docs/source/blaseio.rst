.. module:: blase.blaseio

============================
The write_blender function
============================

The ``write_blender`` function is to run blase and blender in the background.

>>> from ase.build import molecule
>>> from blaseio import write_blender
>>> atoms = molecule('C2H6SO')
>>> batoms = {'label': 'c2h6so', 'atoms': atoms, 'model_type': '1'}
>>> blase = {'output_image': 'figs/c2h6so',}
>>> write_blender(batoms = batoms, blase = blase)

Here, the first keyword ``batoms`` specifies the ASE atomic structure for :class:`Batoms`, and 
second keyword ``blase`` keywords to specify the setting for :class:`Blase`.  Other
possible keywords are: ``display``.

