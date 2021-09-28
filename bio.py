import os
from ase import io
from blase.batoms import Batoms

def read(filename, **kwargs):
    atoms = io.read(filename=filename, **kwargs)
    batoms = Batoms(label = os.path.split(filename)[1], atoms=atoms)
    return batoms