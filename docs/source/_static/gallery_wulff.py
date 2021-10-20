from ase.cluster import wulff_construction
from blase.batoms import Batoms
from blase.bdraw import draw_plane
from blase.butils import removeAll
removeAll()
surfaces = [(1, 1, 1), (1, 0, 0)]
energies = [1.28, 1.69]
atoms = wulff_construction('Au', surfaces, energies, 500, 'fcc')
del atoms[atoms.positions[:, 2] < 0]
nano = Batoms('wulff', atoms = atoms)
nano.show_unit_cell = False
draw_plane(size = 1000, location = (0, 0, min(nano.positions[:, 2]) - nano['Au'].size[0]),  
        color = [0.9, 0.9, 0.9, 1.0])
nano.render.light_type = 'POINT'
nano.render.lock_light_to_camera = False
nano.render.light_energy = 500000
nano.render.light_loc = [30, 4, 30]
nano.render.run([1, -0.3, 0.4], canvas = np.array([[-20, -20, -10], [20, 20, 20]]), 
        engine = 'cycles', output = 'gallery_wulff.png')