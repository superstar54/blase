# from blase.build import molecule
from ase.build import molecule
from ase.visualize import view
ch3choh = molecule('CH3CH2OH')
ch3choh[1].x += 0.4
ch3choh[4].x -= 0.4
ch3choh.write('datas/ch3ch2oh.xyz')
# view(ch3choh)
# from blase.batoms import Batoms
# ch3choh = Batoms('ch3ch2oh', atoms = ch3choh)