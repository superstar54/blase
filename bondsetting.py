"""
"""

class Bondsetting():
    """

    
    """
    def __init__(self, label, bondtalbe = None) -> None:
        self.label = label
    def get_data(self):
        data = {}
        coll = bpy.data.collections[self.label]
        for b in coll.bond:
            data[(b.symbol1, b.symbol2)] = [b.bondlength, b.polyhedra, b.search]
        return data
    @property
    def data(self):
        return self.get_data()
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Bondpair  Bondlength  Polyhedra   Search_bond \n'
        data = self.data
        for key, value in data.items():
            s += '{0:5s} {1:5s} {2:4.3f}       {3:10s}   {4:10s} \n'.format(key[0], key[1], value[0], str(value[1]), str(value[2]))
        s += '-'*60 + '\n'
        return s
    def __getitem__(self, index):
        return self.data[index]
    def __setitem__(self, index, value):
        """
        Add bondpair one by one
        """
        coll = bpy.data.collections[self.label]
        flag = False
        for b in coll.bond:
            if (b.symbol1, b.symbol2) == index:
                b.bondlength = value[0]
                b.polyhedra = value[1]
                b.search = value[2]
                flag = True
        if not flag:
            bond = coll.bond.add()
            bond.symbol1 = index[0]
            bond.symbol2 = index[1]
            bond.bondlength = value[0]
            bond.polyhedra = value[1]
            bond.search = value[2]
