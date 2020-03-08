import os
#
blender_cmd = os.environ['BLENDER_COMMAND']
print(blender_cmd)
errcode = os.system(blender_cmd)
print(errcode)

# blase_path = os.environ['BLASE_PATH']
# blase_cmd = blase_path + '/bin/run-blase.py'
# if display:
#     cmd = blender_cmd + ' -P ' + blase_cmd
# elif queue == 'SLURM':
#     cmd = 'srun -n $SLURM_NTASKS ' +  blender_cmd + ' -b ' + ' -P ' + blase_cmd
# else:
#     cmd = blender_cmd + ' -b ' + ' -P ' + blase_cmd
# # print(cmd)
# errcode = os.system(cmd)