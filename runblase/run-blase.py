from ase.io import read
from ase.io.cube import read_cube_data
from blase.bio import Blase
from blase.batoms import Batoms
import pickle



inputfile = 'blase.inp'

def main():
    print('='*30)
    print('Running blender')
    print('='*30)
    with open(inputfile, 'rb') as f:
        atoms, kwargs = pickle.load(f)
    #
    from blase.bio import Blase
    print('Rendering atoms')
    batoms = Batoms(atoms)
    obj = Blase(batoms, **kwargs)
    obj.render()
    # for function in bobj.functions:
        # name, paras = function
        # getattr(bobj, name)(**paras)
    # if bobj.'run_render':
    # bobj.load_frames()
    
    # bobj.render()
    # bobj.export('h2o.xyz')
    print('\n Finished!')

#
main()

