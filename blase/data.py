material_styles_dict = {
            'jmol'    : {'Specular': 1.0, 'Roughness': 0.001, 'Metallic': 1.0},
            'ase3'    : {'Metallic': 1.0, 'Roughness': 0.001},
            'ceramic' : {'Subsurface': 0.1, 'Metallic': 0.02, 'Specular': 0.5, 'Roughness': 0.0},
            'plastic' : {'Metallic': 0.0, 'Specular': 0.5, 'Roughness': 1.0, },
            'glass'   : {'Metallic': 0.0, 'Specular': 0.5, 'Roughness': 0.0, 'Transmission': 0.98},
            'blase'   : {'Metallic': 0.1, 'Specular': 0.2, 'Roughness': 0.2, },
            'mirror'  : {'Metallic': 0.99, 'Specular': 2.0, 'Roughness': 0.001},
            }
default_batoms = {
        'show_unit_cell': 'default',
        'celllinewidth': 0.025,  # radius of the cylinders representing the cell
        'balltypes': None,
        'radii': None, 
        'colors': None,
        'make_real': False,
        'kind_props': None,
        'bond_cutoff': None,  # 
        'bond_list': {},  # [[atom1, atom2], ... ] pairs of bonding atoms
        'polyhedra_dict': {},
        'search_pbc': False, #{'bonds_dict': {}, 'molecule_list': {}},
        'search_bond': False, #{'bonds_dict': {}, 'molecule_list': {}},
        'search_molecule': False, #{'search_list': None},
        'boundary_list': [],
        'isosurface':None,
        'cube': None,
        'highlight': None, # highlight atoms
}
default_settings = {
        'show_unit_cell': 'default',
        'display': False,  # display while rendering
        'transparent': True,  # transparent background
        'resolution_x': 1000,
        'resolution_y': None,  # 
        'camera': True,
        'camera_loc': None,  # x, y is the image plane, z is *out* of the screen
        'camera_type': 'ORTHO',  #  ['PERSP', 'ORTHO']
        'ortho_scale': None, #
        'camera_lens': 10,  #
        'fstop': 0.5,
        'camera_target': None, #
        'world': False,
        'light': True,
        'light_loc': [0, 0, 200],
        'light_type': 'SUN', # 'POINT', 'SUN', 'SPOT', 'AREA'
        'point_lights': [],  # 
        'light_strength': 1.0,
        'background': 'White',  # color
        'textures': None,  # length of atoms list of texture names
        'engine': 'BLENDER_WORKBENCH', #'BLENDER_EEVEE' #'BLENDER_WORKBENCH'
        'transmits': None,  # transmittance of the atoms
        'bbox': None,
        'bondlinewidth': 0.10,  # radius of the cylinders representing bonds
        'functions': [],
        'run_render': True,
        'animation': False,
        'save_to_blend': False,
        'queue': None,
        'gpu': True,
        'num_samples': 128,
        'build_collection': True,
        }  