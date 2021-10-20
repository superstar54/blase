import pickle
import os

def render(atoms, batoms_input = {}, render_input = {}, display = False, queue = None, ):
    with open('.batoms.inp', 'wb') as f:
        pickle.dump([atoms, batoms_input, render_input], f)
    #
    blender_cmd = 'blender'
    if 'BLENDER_COMMAND' in os.environ.keys():
        blender_cmd = os.environ['BLENDER_COMMAND']
    root = os.path.normpath(os.path.dirname(__file__))
    script = os.path.join(root, 'apirun.py')
    if display:
        cmd = blender_cmd + ' -P ' + script
    elif queue == 'SLURM':
        cmd = 'srun -n $SLURM_NTASKS ' +  blender_cmd + ' -b ' + ' -P ' + script
    else:
        cmd = blender_cmd + ' -b ' + ' -P ' + script
    errcode = os.system(cmd)
    if errcode != 0:
        raise OSError('Command ' + cmd +
                      ' failed with error code %d' % errcode)