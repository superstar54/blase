import pickle
import os
import numpy as np

def write_blender(atoms, batoms = {}, blase = {}, display = False, queue = None, ):
    with open('blase.inp', 'wb') as f:
        pickle.dump([atoms, batoms, blase], f)
    #
    blender_cmd = 'blender'
    if 'BLENDER_COMMAND' in os.environ.keys():
        blender_cmd = os.environ['BLENDER_COMMAND']
    blase_path = os.environ['BLASE_PATH']
    blase_cmd = blase_path + '/runblase/run-blase.py'
    if display:
        cmd = blender_cmd + ' -P ' + blase_cmd
    elif queue == 'SLURM':
        cmd = 'srun -n $SLURM_NTASKS ' +  blender_cmd + ' -b ' + ' -P ' + blase_cmd
    else:
        cmd = blender_cmd + ' -b ' + ' -P ' + blase_cmd
    print(cmd)
    errcode = os.system(cmd)
    # if errcode != 0:
    #     raise OSError('Command ' + cmd +
    #                   ' failed with error code %d' % errcode)




if __name__ == "__main__":
    from ase.io import read, write
    from ase.build import molecule, bulk
    from ase.visualize import view
    from ase.data.colors import jmol_colors
    # test pbc
    # atoms = bulk('Pt', cubic = True)
    atoms = read('../examples/datas/ter-001-la.cif')
    # atoms = atoms*[3, 3, 1]
    # print(atoms)
    # atoms = read('../examples/datas/ceo2.cif')
    # images = []
    # for i in range(10):
        # atoms.positions += np.array([0, 0, 0.2])
        # images.append(atoms.copy())
    #
    # atoms = molecule('CO')
    # atoms.rotate('x', 90)
    # atoms.center(1.0)
    # atoms.pbc = True
    batoms = {
    # 'show_unit_cell': False,
    'model_type': '0',
    # 'scale': 0.5,
    'boundary': [0.05, 0.05, 0.00],
    'remove_bonds': {'La':['O', 'N']},
    # 'polyhedra_dict': {'Ti': ['O', 'N']},
    # 'movie': True,
    'kind_props':{'La': {'color': jmol_colors[26]}},
    }
    blase = {
        #   'camera_loc': camera_loc,  # distance from camera to front atom
        #   'camera_type': 'PERSP',  #  ['PERSP', 'ORTHO', 'PANO']
        #   'camera_lens': 100,  #
        #   'camera_target': atoms.get_center_of_mass(), #
        #   'ortho_scale': None, #
        #   'resolution_x': 1000,
          'output_image': 'lton.png',
        #   'debug': True,
          }
          
    write_blender(atoms, batoms, blase, 
                    display = True,
                )