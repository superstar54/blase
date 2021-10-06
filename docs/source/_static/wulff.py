from ase.cluster import wulff_construction
from blase.batoms import Batoms
from blase.bdraw import draw_plane

surfaces = [(1, 1, 1), (1, 0, 0)]
energies = [1.28, 1.69]
atoms = wulff_construction('Au', surfaces, energies, 5000, 'fcc')
atoms.center(vacuum=2.0)
# view(atoms)
camera_loc =  atoms.get_center_of_mass() + [0, -300, 200]

nano = Batoms('wulff', atoms = atoms)
nano.show_unit_cell = False
draw_plane(size = 1000, location = (0, 0, -1.0),  color = [1.0, 1.0, 1.0, 1.0], material_style = 'mirror')
nano.render.light_energy = 25
nano.render.run(engine = 'cycles', camera_loc = camera_loc, camera_type = 'PERSP', camera_lens = 100, 
                camera_target = atoms.get_center_of_mass(), num_samples = 32,
                output = 'wulff.png'
                )