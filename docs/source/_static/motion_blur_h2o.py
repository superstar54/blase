from blase.batoms import Batoms
from ase import Atoms
from ase.build import molecule

h2o = molecule('H2O')
images = []
for i in range(20):
    atoms = Atoms()
    for j in [-1, 1]:
        for k in [-1, 1]:
            temp = h2o.copy()
            temp.translate([i/3*j, i*i/50*k, 0])
            atoms = atoms + temp
    images.append(atoms)
h2o = Batoms(label = 'h2o', atoms = images)
h2o.load_frames()
h2o.render.run(engine = 'eevee', frame = 10, use_motion_blur = True, ortho_scale = 30, ratio = 1, transparent = False)