import bpy
import addon_utils
import sys
import os
import subprocess
#bpy.ops.wm.addon_install(filepath='/home/shane/Downloads/testaddon.py')
#path = os.path.join(sys.prefix, 'bin', 'python3.9')
#subprocess.call([path, "-m", "pip", "install", "ase"])
#subprocess.call([path, "-m", "pip", "install", "scikit-image"])


for addon in bpy.context.preferences.addons:
    print(addon.module)

addon_utils.enable('blase', default_set=True) 

for addon in bpy.context.preferences.addons:
    print(addon.module)
