#!/usr/bin/env python
from ase.io import read
from ase.io.cube import read_cube_data
from blase.bio import Blase
import pickle



inputfile = 'blase.inp'

def main():
    print('='*30)
    print('Running blender')
    print('='*30)
    with open(inputfile, 'rb') as f:
        images, kwargs = pickle.load(f)
    #
    bobj = Blase(images, **kwargs)
    # view(images)
    bobj.draw_cell()
    bobj.draw_atoms()
    bobj.draw_bonds()
    for function in bobj.functions:
        name, paras = function
        getattr(bobj, name)(**paras)
    # if bobj.'run_render':
    # bobj.load_frames()
    bobj.render()
    print('\n Finished!')

#
main()

