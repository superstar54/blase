from blase.atoms import Batoms
atoms = Batoms('H2O', positions = [[0, 0 ,0], [0, 0, 1], [0, 0, 3]], name = 'h2o')
atoms.draw()
atoms.delete([1, 2])
