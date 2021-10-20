#!/usr/bin/env python
from ase.io import read
from ase.build import molecule
from ase.io.cube import read_cube_data
import pickle
import sys
import argparse
import os


def save(batoms_input, render_input):
    with open('.batoms.inp', 'wb') as f:
        pickle.dump([batoms_input, render_input], f)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputfile', type=str, default='',
                        help="the input json file, includes coordinates of a \
                        set of points, threshold and a list of pairs of points ")
    parser.add_argument('--display', action='store_true', default=True,
                        help="render")
    parser.add_argument('--run_render', action='store_false', default=False,
                        help="render")
    parser.add_argument('--model_type', '-m', type=str, default='0',
                        help="structure model")
    parser.add_argument('--level', '-iso', type=str, default='0.002',
                        help="structure model")
    parser.add_argument('--output_image', '-o', type=str, default='output',
                        help="write output to specified file ")
    parser.add_argument('--camera', action='store_false', default=False,
                        help="camera")
    parser.add_argument('--light', action='store_false', default=False,
                        help="light")
    args = parser.parse_args()
    #
    render_input = {}
    batoms_input = {}
    render_input['output_image'] = args.output_image
    render_input['run_render'] = args.run_render
    batoms_input['inputfile'] = args.inputfile
    batoms_input['model_type'] = args.model_type
    save(batoms_input, render_input)
    #-----------
    root = os.path.normpath(os.path.dirname(__file__))
    script = os.path.join(root, 'run.py')
    #
    blender_cmd = 'blender'
    if 'BLENDER_COMMAND' in os.environ.keys():
        blender_cmd = os.environ['BLENDER_COMMAND']
    if args.display:
        cmd = blender_cmd + ' -P ' + script
    else:
        cmd = blender_cmd + ' -b ' + ' -P ' + script
    errcode = os.system(cmd)

if __name__ == "__main__":
    main()