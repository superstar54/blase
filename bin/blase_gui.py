#!/usr/bin/env python
from ase.io import read
from ase.io.cube import read_cube_data
from blase.bio import Blase
from blase.btools import draw_cell, draw_atoms, draw_bonds, draw_polyhedras, draw_isosurface, bond_source, cylinder_mesh_from_instance, clean_default
import sys


def import_blase(inputfile, 
               model_type = '0',
               camera = 'True',
               light = 'True',
               world = 'False',
               show_unit_cell = False,
               ball_type = 'Ball-and-stick',
               radii = 1.0,
               bond_cutoff = 1.0,
               search_pbc = False,
               search_molecule = False,
               make_real = False,
               ):
    #
    print('='*30)
    print('Import structure')
    print('='*30)
    print(model_type)
    #
    kwargs = {'show_unit_cell': show_unit_cell, 
          'radii': radii,
          'bond_cutoff': bond_cutoff,
          'world': world,
          'camera': camera,
          'light': light,
          'search_pbc':search_pbc,
          'search_molecule': search_molecule,
          'make_real': make_real,
          'display': True,
          }
    
    if inputfile.split('.')[-1] == 'cube':
        images, data = read_cube_data(inputfile)
    else:
        images = read(inputfile)
    
    if model_type == '0':
        # view(images)
        kwargs['radii'] = 0.6
        kwargs['bond_cutoff'] = 1.0
        bobj = Blase(images, **kwargs)
        draw_cell(bobj)
        draw_atoms(bobj)
        draw_bonds(bobj)
    elif model_type == '1':
        kwargs['bond_cutoff'] = None
        bobj = Blase(images, **kwargs)
        draw_cell(bobj)
        draw_atoms(bobj)
    elif model_type == '2':
        kwargs['radii'] = 0.6
        kwargs['bond_cutoff'] = 1.0
        kwargs['polyhedral'] = {}
        bobj = Blase(images, **kwargs)
        draw_cell(bobj)
        draw_atoms(bobj)
        draw_bonds(bobj)
        draw_polyhedras(bobj)
    elif model_type == '3':
        # kwargs['bond_cutoff'] = 1.0
        bobj = Blase(images, **kwargs)
        draw_cell(bobj)
        draw_bonds(bobj)

    print('\n Finished!')

#
inputfile = sys.argv[3]
model_type = sys.argv[4]

import_blase(inputfile, model_type)

