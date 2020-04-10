#!/usr/bin/env python
from ase.io import read
from ase.io.cube import read_cube_data
from blase.bio import Blase
from blase.btools import draw_cell, draw_atoms, draw_bonds, draw_polyhedras, draw_isosurface, bond_source, cylinder_mesh_from_instance, clean_default
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
    bobj.draw()
    for function in bobj.functions:
        name, paras = function
        getattr(bobj, name)(**paras)
    # if bobj.'run_render':
    # bobj.load_frames()
    
    bobj.render()
    # bobj.export('h2o.xyz')
    print('\n Finished!')

#
main()

