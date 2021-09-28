from ase.io import read
from ase.io.cube import read_cube_data
from blase.render import Blase
from blase.batoms import Batoms
from blase.tools import get_bbox
import pickle



inputfile = 'blase.inp'

def main():
    with open(inputfile, 'rb') as f:
        batoms, blase = pickle.load(f)
    #
    if isinstance(batoms, dict):
        batoms = [batoms]
    bboxs = []
    for batom in batoms:
        ba = Batoms(**batom)
        bbox = get_bbox(bbox = None, atoms = batom['atoms'])
        bboxs.append(bbox)
        ba.draw()
        # batoms.load_frames()
    print('--------------Render--------------')
    print('Rendering atoms')
    if 'bbox' not in blase:
        bbox = bboxs[0]
        for xyz in bboxs[1:]:
            for i in range(3):
                bbox[i][0] = min(bbox[i][0], xyz[i][0])
                bbox[i][1] = max(bbox[i][1], xyz[i][1])
        blase['bbox'] = bbox
    # print(blase)
    bobj = Blase(**blase)
    for function in bobj.functions:
        name, paras = function
        getattr(bobj, name)(**paras)
    bobj.render()
    print('-'*20)
    print('\n Finished!')

#
main()

