from ase.build import molecule
from blase.batoms import Batoms
ch3choh = molecule('CH3CH2OH')
ch3choh = Batoms('ch3ch2oh', atoms = ch3choh)
ch3choh.bondsetting['H-O'].min = 2.0
ch3choh.bondsetting['H-O'].max = 3.0
ch3choh.bondsetting['H-O'].width = 0.01
ch3choh.bondsetting['H-O'].style = '1'
ch3choh.model_type = 1

