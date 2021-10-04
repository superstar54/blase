from blase.batoms import Batoms
a = 3.96
positions = [[0, 0, 0], [a/2, a/2, 0], [a/2, 0, a/2], [0, a/2, a/2]]
pt = Batoms({'Pt': positions}, pbc = True, cell = (a, a, a))
pt.segments = 6
pt.repeat([50, 50, 50])
pt.render()
