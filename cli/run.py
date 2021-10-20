from ase.io import read
from ase.io.cube import read_cube_data
from blase.batoms import Batoms
import pickle
import os

def main():
    with open('.batoms.inp', 'rb') as f:
        batoms_input, render_input = pickle.load(f)
        inputfile = batoms_input['inputfile']
        head, tail = os.path.split(os.path.split(inputfile)[0])
        label = tail.split('.')[0]
        if tail.split('.')[-1] == 'cube':
            volume, atoms = read_cube_data(inputfile)
            batoms = Batoms(label, atoms = atoms, volume=volume)
        else:
            atoms = read(inputfile)
            batoms = Batoms(label = label, atoms = atoms)
    batoms.model_type = batoms_input['model_type']
    if atoms.pbc.any():
      batoms_input['show_unit_cell'] = True
    batoms.render.run(output_image = render_input['output_image'], 
            run_render = render_input['run_render'])

if __name__ == "__main__":
    main()
