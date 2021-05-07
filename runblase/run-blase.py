from ase.io import read
from ase.io.cube import read_cube_data
from blase.bio import Blase
from blase.batoms import Batoms
import pickle



inputfile = 'blase.inp'

def main():
    with open(inputfile, 'rb') as f:
        atoms, batoms, blase = pickle.load(f)
    #
    batoms = Batoms(atoms, **batoms)
    batoms.draw()
    # batoms.load_frames()
    batoms.render(**blase)
    print('-'*20)
    print('\n Finished!')

#
main()

