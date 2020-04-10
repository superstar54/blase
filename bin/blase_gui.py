#!/usr/bin/env python
from ase.io import read
from ase.io.cube import read_cube_data
from blase.tools import write_blender, get_bondpairs
import sys
import argparse

# default
kwargs = {
      'radii': 1.0,
      'bond_cutoff': 1.0,
      }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputfile', type=str,
                        help="the input json file, includes coordinates of a \
                        set of points, threshold and a list of pairs of points ")
    parser.add_argument('--display', action='store_true', default=True,
                        help="render")
    parser.add_argument('--run_render', action='store_false', default=False,
                        help="render")
    parser.add_argument('--model', '-m', type=int, default=0,
                        help="structure model")
    parser.add_argument('--outfile', '-o', type=str, default='output',
                        help="write output to specified file ")
    parser.add_argument('--camera', action='store_false', default=False,
                        help="camera")
    parser.add_argument('--light', action='store_false', default=False,
                        help="light")
    parser.add_argument('--search_pbc_atoms', action='store_false', default=False,
                        help="search_pbc_atoms")
    parser.add_argument('--search_molecule', action='store_false', default=False,
                        help="search_molecule")
    args = parser.parse_args()
    #
    if args.inputfile.split('.')[-1] == 'cube':
        atoms, data = read_cube_data(inputfile)
    else:
        atoms = read(args.inputfile)
    if atoms.pbc.any():
      kwargs['show_unit_cell'] = True
    kwargs['display'] = args.display
    kwargs['camera'] = args.camera
    kwargs['light'] = args.light
    kwargs['outfile'] = args.outfile+'.png'
    kwargs['run_render'] = args.run_render
    kwargs['search_pbc_atoms'] = args.search_pbc_atoms
    kwargs['search_molecule'] = args.search_molecule
    
    if args.model == 0:
        kwargs['radii'] = 0.6
        kwargs['bond_cutoff'] = 1.0
    elif args.model == 1:
        kwargs['bond_cutoff'] = None
    elif args.model == 2:
        kwargs['radii'] = 0.6
        kwargs['bond_cutoff'] = 1.0
        kwargs['polyhedral'] = {}
    
    print('='*30)
    print('Import structure')
    print('='*30)
    write_blender(atoms, **kwargs)
    print('\n Finished!')


if __name__ == "__main__":
    main()