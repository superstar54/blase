import bpy

def object_mode():
    for object in bpy.data.objects:
            if object.mode == 'EDIT':
                bpy.ops.object.mode_set(mode = 'OBJECT')