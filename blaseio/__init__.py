import pickle
import os
import numpy as np

def write_blender(batoms = [], blase = {}, display = False, queue = None, ):
    with open('blase.inp', 'wb') as f:
        pickle.dump([batoms, blase], f)
    #
    blender_cmd = 'blender'
    if 'BLENDER_COMMAND' in os.environ.keys():
        blender_cmd = os.environ['BLENDER_COMMAND']
    blase_path = os.environ['BLASE_PATH']
    blase_cmd = os.path.join(blase_path, 'run-blase.py')
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
