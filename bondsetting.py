"""
"""
import bpy
from blase.btools import object_mode

class Bondsetting():
    """
    Set bond infomation.
    
    """
    def __init__(self, label, bondtalbe = None) -> None:
        self.label = label
    def get_data(self):
        data = {}
        coll = bpy.data.collections[self.label]
        for b in coll.bond:
            data[(b.symbol1, b.symbol2)] = [b.min, b.max, b.polyhedra, b.search]
        return data
    @property
    def data(self):
        return self.get_data()
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Bondpair      min     max   Polyhedra   Search_bond \n'
        data = self.data
        for key, value in data.items():
            s += '{0:5s} {1:5s} {2:4.3f}   {3:4.3f}      {4:10s}   {5:10s} \n'.format(key[0], key[1], value[0], value[1], str(value[2]), str(value[3]))
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
                b.min = value[0]
                b.max = value[1]
                b.polyhedra = value[2]
                b.search = value[3]
                flag = True
        if not flag:
            bond = coll.bond.add()
            bond.symbol1 = index[0]
            bond.symbol2 = index[1]
            bond.min = value[0]
            bond.max = value[1]
            bond.polyhedra = value[2]
            bond.search = value[3]
    def copy(self, label):
        object_mode()
        bondsetting = Bondsetting(label)
        for key, value in self.data.items():
            bondsetting[key] = value
        return bondsetting