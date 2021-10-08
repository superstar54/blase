.. module:: blase.polyhedrasetting

===========================
The Polyhedrasetting object
===========================

The :class:`Polyhedrasetting` object is used to store and set all parameters related with polyhedra. It is a collection of :class:`BlasePolyhedra` object. It should always bind with a :class:`Batoms` object. Possible keywords are: ``symbol``, ``color`` and ``linewidth``. 



You can print the default polyhedrasetting by:

>>> ch4.polyhedrasetting

One can change color for polyhedra ``C`` by. 

>>> ch4.polyhedrasetting['C'].color = [0.8, 0.1, 0.3, 0.3]
>>> ch4.model_type = 2




List of all Methods
===================

.. autoclass:: Polyhedrasetting
   :members: