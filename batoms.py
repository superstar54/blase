"""
Definition of the Batoms class in the blase package.

"""

from numpy.core.fromnumeric import shape
import bpy
from ase import Atoms
from blase.batom import Batom
from blase.bondsetting import Bondsetting, build_bondlists, search_skin, build_polyhedralists, \
                              get_bondtable
from blase.cell import Bcell
from blase.render import Render   
from blase.boundary import search_boundary
from blase.bdraw import draw_cell_cylinder, draw_bond_kind, draw_polyhedra_kind, \
                        draw_isosurface
from blase.butils import object_mode
import numpy as np
from time import time

import logging
logging.basicConfig(
                    format=('%(levelname)-8s '
                            '[%(funcName)-20s]: %(message)s'),
                    level=logging.INFO)

logger = logging.getLogger(__name__)

subcollections = ['atom', 'cell', 'bond', 'polyhedra', 'instancer', 
'instancer_atom', 'isosurface', 'virtual', 'boundary', 'skin', 'render', 'text']

class Batoms():
    """
    Batoms object

    The Batoms object is a interface to a blase collection in Blender.
    A blase collections is organised in the following way, 
    take water molecule as a example:

    * h2o                          # main collection

      * h2o_atom                   # atoms collection
    
        * atom_h2o_H               # atoms object
    
        * atom_h2o_O               # atoms object

      * h2o_bond                   # bonds collection

        * bond_h2o_H               # bond object
    
        * bond_h2o_O               # bond object
    
      * h2o_cell                   # cell collection
    
      * h2o_instancer              # instancer collection
    
        * instancer_atom_h2o_H     # sphere object
    
        * instancer_atom_h2o_O     # sphere object
    
    Then, a Batoms object is linked to this main collection in Blender. 

    Parameters:

    label: str
        Name for the collection in Blender.
    species: dict or list
        Can be a dict with symbols and positions. Examples:
        {
         'O': [[0, 0, 0.40]], 
         'H': [[0, -0.76, -0.2], [0, 0.76, -0.2]]
        }
        Or can be a list of Baom object.
        [Batom('h2o', 'H', ...), Batom('h2o', 'O', ...)]
    atoms: ase.atoms.Atoms object or a list of ase.atoms.Atoms object
        or pymatgen structure object
    model_type: str
        enum in '0', '1', '2', '3', Default value: '0'
    pbc: Bool or three Bool
        Periodic boundary conditions. Examples: True,
        False, (True, False, False).  Default value: False.
    cell: 3x3 matrix or length 3 or 6 vector
        Unit cell.
    segments: list
        value should be int, and in [3, 100000]
        segments and ring_count in bpy.ops.mesh.primitive_uv_sphere_add

    boundary:  list 
        search atoms at the boundary

    Examples:
    >>> from blase import Batoms
    >>> h2o = Batoms(label = 'h2o', species = {'O': [[0, 0, 0.40]], 
    ...             'H': [[0, -0.76, -0.2], [0, 0.76, -0.2]]})
    
    Here is equivalent:
    
    >>> h = Batom(label = 'h2o', species = 'H', 
    ...           positions = [[0, -0.76, -0.2], [0, 0.76, -0.2]])
    >>> o = Batom(label = 'h2o', species = 'O', 
    ...           positions = [[0, 0, 0.40]])
    >>> h2o = Batoms('h2o', [h, o])

    """
    

    def __init__(self, label = None,
                species = None,
                atoms = None, 
                pbc = False, cell = None,
                bondsetting = None,
                render = None,
                model_type = '0', 
                boundary = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
                show_unit_cell = True,
                isosurface = [],
                segments = [32, 16],
                shape = 'UV_SPHERE',
                kind_props = {},
                color_style = 'JMOL',
                material_style = 'blase',
                bsdf_inputs = None,
                movie = False,
                draw = True, 
                 ):
        #
        self.scene = bpy.context.scene
        self.bondsetting = bondsetting
        self.render = render
        self.segments = segments
        self.shape = shape
        self.isosurface = isosurface
        self.kind_props = kind_props
        self.label = label
        self.color_style = color_style
        self.material_style = material_style
        self.bsdf_inputs = bsdf_inputs
        if species:
            if not self.label:
                self.label = ''.join(['%s%s'%(species, len(positions)) for sp, positions in species.items()])
            self.set_collection(model_type, boundary)
            self.from_species(species, pbc, cell)
        elif atoms:
            if isinstance(atoms, list):
                self.images = atoms
                atoms = self.images[0]
            self.set_collection(model_type, boundary)
            if 'ase' in str(type(atoms)):
                self.from_ase(atoms)
            elif 'pymatgen' in str(type(atoms)):
                self.from_pymatgen(atoms)
        elif self.label:
            print('Build from collection')
            self.from_collection(self.label)
        else:
            raise Exception("Failed, species, atoms or coll should be provided!"%self.label)
        if not self.bondsetting:
            self.bondsetting = Bondsetting(self.label)
            bondtable = get_bondtable(self.species, cutoff=1.2)
            for key, value in bondtable.items():
                self.bondsetting[key] = value
        self.coll.blase.show_unit_cell = show_unit_cell
        if not self.render:
            self.render = Render(self.label, batoms = self)
        if draw:
            self.draw()
        if movie:
            self.load_frames()
        self.show_index()
    def from_species(self, species, pbc = None, cell = None):
        """
        """
        if isinstance(species, dict):
            for sp, positions in species.items():
                if sp not in self.kind_props: self.kind_props[sp] = {}
                ba = Batom(self.label, sp, positions, segments = self.segments, shape = self.shape,
                            props = self.kind_props[sp], material_style=self.material_style, 
                            bsdf_inputs=self.bsdf_inputs, color_style=self.color_style)
                self.coll.children['%s_atom'%self.label].objects.link(ba.batom)
                self.coll.children['%s_instancer'%self.label].objects.link(ba.instancer)
        elif isinstance(species, list):
            for batom in species:
                if not isinstance(batom, Batom):
                    raise Exception('%s is not a Batom object.'%batom)
                self.coll.children['%s_atom'%self.label].objects.link(batom.batom)
                self.coll.children['%s_instancer'%self.label].objects.link(batom.instancer)
        self.coll.is_batoms = True
        self._cell = Bcell(self.label, cell)
        self.coll.children['%s_cell'%self.label].objects.link(self._cell.bcell)
        self.set_pbc(pbc)
    def from_ase(self, atoms):
        """
        Import structure from ASE atoms.
        """
        if 'species' not in atoms.info:
            atoms.info['species'] = atoms.get_chemical_symbols()
        species_list = list(set(atoms.info['species']))
        for species in species_list:
            indices = [index for index, x in enumerate(atoms.info['species']) if x == species]
            if species not in self.kind_props: self.kind_props[species] = {}
            ba = Batom(self.label, species, atoms.positions[indices], 
                        segments = self.segments, shape = self.shape, props = self.kind_props[species], 
                        material_style=self.material_style, bsdf_inputs=self.bsdf_inputs, 
                        color_style=self.color_style)
            self.coll.children['%s_atom'%self.label].objects.link(ba.batom)
            self.coll.children['%s_instancer'%self.label].objects.link(ba.instancer)
        self.coll.is_batoms = True
        self.coll.blase.pbc = self.npbool2bool(atoms.pbc)
        self._cell = Bcell(self.label, atoms.cell)
        self.coll.children['%s_cell'%self.label].objects.link(self._cell.bcell)
    def from_pymatgen(self, structure):
        """
        Import structure from Pymatgen structure.
        """
        symbols = [str(site.specie.symbol) for site in structure]
        if hasattr(structure, "lattice"):
            cell = structure.lattice.matrix
            pbc = True
        else:
            cell = None
            pbc = False
        species_list = list(set(symbols))
        for species in species_list:
            positions = [structure[index].coords for index, x in enumerate(symbols) if x == species]
            ba = Batom(self.label, species, positions, segments = self.segments, shape = self.shape, material_style=self.material_style, bsdf_inputs=self.bsdf_inputs, color_style=self.color_style)
            self.coll.children['%s_atom'%self.label].objects.link(ba.batom)
            self.coll.children['%s_instancer'%self.label].objects.link(ba.instancer)
        self.coll.is_batoms = True
        self.set_pbc(pbc)
        self._cell = Bcell(self.label, cell)
        self.coll.children['%s_cell'%self.label].objects.link(self._cell.bcell)
    def from_collection(self, collection_name):
        """
        """
        if collection_name not in bpy.data.collections:
            raise Exception("%s is not a collection!"%collection_name)
        elif not bpy.data.collections[collection_name].is_batoms:
            raise Exception("%s is not Batoms collection!"%collection_name)
        self.label = collection_name
        self._cell = Bcell(label = collection_name)
        self.bondsetting = Bondsetting(self.label)

    def npbool2bool(self, pbc):
        """
        """
        newpbc = []
        for i in range(3):
            if pbc[i]:
                newpbc.append(True)
            else:
                newpbc.append(False)
        return newpbc
    def set_collection(self, model_type = '0', boundary = [0, 0, 0]):
        """
        build main collection and its child collections.
        """
        for coll in bpy.data.collections:
            if self.label == coll.name:
                raise Exception("Failed, the name %s already in use!"%self.label)
        coll = bpy.data.collections.new(self.label)
        self.coll.blase.is_batoms = True
        self.scene.collection.children.link(self.coll)
        for sub_name in subcollections:
            subcoll = bpy.data.collections.new('%s_%s'%(self.label, sub_name))
            self.coll.children.link(subcoll)
        self.coll.blase.model_type = model_type
        self.coll.blase.boundary = boundary
    
    def draw_cell(self, celllinewidth = 0.03):
        """
        Draw unit cell
        """
        object_mode()
        cell_vertices = self.cell.verts
        if np.max(abs(cell_vertices)) < 1e-6:
            return 0
        self.clean_blase_objects('cell', ['cylinder', 'point'])
        if self.show_unit_cell:
            draw_cell_cylinder(self.coll.children['%s_cell'%self.label], 
                            cell_vertices, 
                            label = self.label, 
                            celllinewidth = celllinewidth)
    def draw_bonds(self, bondlinewidth = 0.1):
        """
        Draw bonds.

        Parameters:

        cutoff: float
            cutoff used to build bond pairs.
        """
        # if not self.bondlist:
        object_mode()
        atoms, n1, n2, n3 = self.get_atoms_with_boundary()
        self.bondlist = build_bondlists(atoms, self.bondsetting.data)
        self.calc_bond_data(atoms, self.bondlist)
        for species, bond_data in self.bond_kinds.items():
            draw_bond_kind(species, bond_data, label = self.label, 
                        coll = self.coll.children['%s_bond'%self.label], bondlinewidth = bondlinewidth)
    def draw_polyhedras(self):
        """
        Draw bonds.
        Parameters:

        cutoff: float
            cutoff used to build bond pairs.
        """
        object_mode()
        atoms, n1, n1, n3 = self.get_atoms_with_boundary()
        polyhedra_kinds = build_polyhedralists(atoms, self.bondlist, self.bondsetting.data, color_style = self.color_style)
        for species, polyhedra_data in polyhedra_kinds.items():
            draw_polyhedra_kind(species, polyhedra_data, label = self.label,
                        coll = self.coll.children['%s_polyhedra'%self.label])
    def draw_isosurface(self, isosurface = []):
        """
        Draw bonds.

        Parameters:

        isosurface: list
            isosurface data.
        """
        object_mode()
        if not isosurface:
            isosurface = self.isosurface
        volume = self.isosurface[0]
        icolor = 0
        if len(self.isosurface) == 1:
            draw_isosurface(self.coll.children['%s_isosurface'%self.label], volume, cell = self.atoms.cell, level=None, icolor = icolor)
        for level in self.isosurface[1:]:
            draw_isosurface(self.coll.children['%s_isosurface'%self.label], volume, cell = self.atoms.cell, level=level, icolor = icolor)
            icolor += 1
    def draw_cavity(self, radius):
        """
        cavity
        for porous materials
        >>> from ase.io import read
        >>> atoms = read('docs/source/_static/datas/mof-5.cif')
        """
        from blase.tools import find_cage
        object_mode()
        self.clean_blase_objects('virtual')
        positions = find_cage(self.cell, self.atoms.positions, radius)
        ba = Batom(self.label, 'Au_cavity', positions, scale = radius/2.8, material_style='blase', bsdf_inputs=self.bsdf_inputs, color_style=self.color_style)
        self.coll.children['%s_virtual'%self.label].objects.link(ba.batom)
        self.coll.children['%s_virtual'%self.label].objects.link(ba.instancer)
    def clean_blase_objects(self, coll, objs = None):
        """
        remove all bond object in the bond collection
        """
        if not objs:
            for obj in self.coll.children['%s_%s'%(self.label, coll)].all_objects:
                bpy.data.objects.remove(obj, do_unlink = True)
        else:
            for obj in objs:
                name = '%s_%s_%s'%(coll, self.label, obj)
                if name in self.coll.children['%s_%s'%(self.label, coll)].all_objects:
                    obj = self.coll.children['%s_%s'%(self.label, coll)].all_objects[name]
                    bpy.data.objects.remove(obj, do_unlink = True)
                
    def show_index(self, index_type = 0):
        """
        """
        bpy.context.preferences.view.show_developer_ui = True
        for a in bpy.context.screen.areas:
            if a.type == 'VIEW_3D':
                overlay = a.spaces.active.overlay
                overlay.show_extra_indices = True
    @property
    def scale(self):
        return self.get_scale()
    @scale.setter
    def scale(self, scale):
        self.set_scale(scale)
    def get_scale(self):
        scale = {}
        for coll in [self.batoms, self.batoms_boundary, self.batoms_skin]:
            for batom in coll.values():
                scale[batom.species] = batom.scale
        return scale
    def set_scale(self, scale):
        for coll in [self.batoms, self.batoms_boundary, self.batoms_skin]:
            for batom in coll.values():
                batom.scale = scale
    def draw(self, model_type = None):
        """
        Draw atoms, bonds, polyhedra.

        Parameters:

        model_type: str

        """
        if model_type is not None and model_type not in ['0', '1', '2', '3']:
            raise Exception('model_type %s should be: 0, 1, 2, 3'%model_type)
        if not model_type:
            model_type = self.model_type
        else:
            self.model_type = model_type
        # self.draw_cell()
        bpy.ops.ed.undo_push()
        self.clean_blase_objects('bond')
        bpy.ops.ed.undo_push()
        self.clean_blase_objects('polyhedra')
        bpy.ops.ed.undo_push()
        if model_type == '0':
            self.scale = 1.0
        elif model_type == '1':
            self.scale = 0.5
            self.draw_bonds()
        elif model_type == '2':
            self.scale = 0.5
            self.draw_bonds()
            self.draw_polyhedras()
        elif model_type == '3':
            self.scale = 0.01
            self.draw_bonds(bondlinewidth = 0.02)
        if self.isosurface:
            self.draw_isosurface()
    def replace(self, species1, species2, index = []):
        """
        replace atoms.

        Parameters:
        
        index: list
            index of atoms will be replaced.

        species1: str

        species2: str
            atoms will be changed to this element.

        >>> from ase.build import molecule, fcc111
        >>> from blase.batoms import Batoms
        >>> pt111 = fcc111('Pt', (5, 5, 4), vacuum = 5.0)
        >>> pt111 = Batoms(atoms = pt111, label = 'pt111')
        >>> pt111.replace('Pt', 'Au', [93])
        >>> pt111.replace('Pt', 'Au', range(20))

        """
        # if kind exists, merger, otherwise build a new kind and add.
        object_mode()
        name = 'atom_%s_%s'%(self.label, species2)
        positions = [self.batoms[species1][i] for i in index]
        if species2 in self.batoms:
            self.batoms[species2].add_vertices(positions)
        else:
            ba = Batom(self.label, species2, positions, segments = self.segments, shape = self.shape, material_style=self.material_style, bsdf_inputs=self.bsdf_inputs, color_style=self.color_style)
            self.coll.children['%s_atom'%self.label].objects.link(ba.batom)
            self.coll.children['%s_instancer'%self.label].objects.link(ba.instancer)
        self.batoms[species1].delete(index)
            
    
    def delete(self, species, index = []):
        """
        delete atoms.

        species: str

        index: list
            index of atoms to be delete
        
        For example, delete the second atom in H species.
        Please note that index start from 0.

        >>> h2o.delete([1])

        """
        self.batoms[species].delete(index)
    def translate(self, displacement):
        """Translate atomic positions.

        The displacement argument is an xyz vector.

        For example, move h2o molecule by a vector [0, 0, 5]

        >>> h2o.translate([0, 0, 5])

        """
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        for obj in self.coll.all_objects:
            if 'instancer' not in obj.name and 'boundary' not in obj.name:
                obj.select_set(True)
        bpy.ops.transform.translate(value=displacement)
    def rotate(self, angle, axis = 'Z', orient_type = 'GLOBAL'):
        """Rotate atomic based on a axis and an angle.

        Parameters:

        angle: float
            Angle that the atoms is rotated around the axis.
        axis: str
            'X', 'Y' or 'Z'.

        For example, rotate h2o molecule 90 degree around 'Z' axis:
        
        >>> h2o.rotate(90, 'Z')

        """
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        for coll in self.coll.children:
            for obj in coll.objects:
                obj.select_set(True)
        bpy.ops.transform.rotate(value=angle, orient_axis=axis.upper(), orient_type = orient_type)
    
    def __getitem__(self, species):
        """Return a subset of the Batom.

        species -- str, describing which batom to return.
        """

        if isinstance(species, str):
            if species not in self.batoms:
                raise SystemExit('%s is not in this structure'%species)
            return self.batoms[species]
        elif isinstance(species, list):
            raise SystemExit('dict not supported yet!')
    def __len__(self):
        return len(self.positions)
    def __add__(self, other):
        self += other
        return self
    def __iadd__(self, other):
        self.extend(other)
        return self
    def __repr__(self) -> str:
        text = []
        text.append('label={0}, '.format(self.label))
        text.append('species='.format(self.cell))
        text.append('%s '%(list(self.batoms)))
        text.append('cell={0}, '.format(self.cell))
        text.append('pbc={0}'.format(self.pbc))
        text = "".join(text)
        text = "Batoms(%s)"%text
        return text
    def extend(self, other):
        """
        Extend batoms object by appending batoms from *other*.
        
        >>> from ase.build import molecule, fcc111
        >>> from blase.batoms import Batoms
        >>> import numpy as np
        >>> co = molecule('CO')
        >>> co = Batoms(atoms = co, draw = True)
        >>> au = fcc111('Au', (5, 5, 4), vacuum=5.0)
        >>> au = Batoms(atoms = au, draw = True)
        >>> co.translate(au.atoms[-1].position + np.array([0, 0, 2]))
        >>> au.extend(co)
        >>> au.write('au111-co.cif')
        
        or,

        >>> au = au + co

        """
        for species, batom in other.batoms.items():
            if species in self.species:
                self.batoms[species].extend(batom)
            else:
                ba = batom.copy(self.label, species)
                t = self.cell.location - ba.location
                ba.batom.location = self.cell.location
                bpy.context.view_layer.update()
                ba.positions = ba.positions - t
                self.coll.children['%s_atom'%self.label].objects.link(ba.batom)
        for key, value in other.bondsetting.data.items():
            self.bondsetting[key] = value
        self.remove_collection(other.label)

    def remove_collection(self, name):
        collection = bpy.data.collections.get(name)
        for obj in collection.all_objects:
            bpy.data.objects.remove(obj, do_unlink=True)
        for coll in collection.children:
            bpy.data.collections.remove(coll)
        bpy.data.collections.remove(collection)
    def __imul__(self, m):
        """
        """
        for species, batom in self.batoms.items():
            batom.repeat(m, self.cell)
        self.cell.repeat(m)
        return self
    def repeat(self, rep):
        """
        Create new repeated atoms object.

        >>> from ase.build import bulk
        >>> from blase.batoms import Batoms
        >>> au = bulk('Au', cubic = True)
        >>> au = Batoms(atoms = au)
        >>> au.draw()
        >>> au.repeat([2, 2, 2])
        
        """
        self.__imul__(rep)
    def __mul__(self, rep):
        self.repeat(rep)
        return self
    def copy(self, label = None):
        """
        Return a copy.

        name: str
            The name of the copy.

        For example, copy h2o molecule:
        
        >>> h2o_new = h2o.copy(label = 'h2o_new')

        """
        if not label:
            label = self.label + 'copy'
        species_dict = {x:self.batoms[x].copy(label, x) for x in self.species}
        batoms = self.__class__(species_dict = species_dict, label = label, cell = self.cell.verts, pbc = self.pbc, model_type = self.coll.blase.model_type)
        batoms.translate([2, 2, 2])
        batoms.bondsetting = self.bondsetting.copy(label)
        return batoms
    def write(self, filename, local = True):
        """
        Save atoms to file.

        >>> h2o.write('h2o.cif')
        
        """
        atoms = self.batoms2atoms(self.batoms, local = local)
        atoms.write(filename)
    def update(self):
        """
        """
        pass
    @property
    def coll(self):
        return self.get_coll()
    def get_coll(self):
        return bpy.data.collections[self.label]
    @property
    def coll_atom(self):
        return self.get_coll_atom()
    def get_coll_atom(self):
        return self.coll.children['%s_atom'%self.label]
    @property
    def coll_boundary(self):
        return self.get_coll_boundary()
    def get_coll_boundary(self):
        return self.coll.children['%s_boundary'%self.label]
    @property
    def coll_skin(self):
        return self.get_coll_skin()
    def get_coll_skin(self):
        return self.coll.children['%s_skin'%self.label]
    @property
    def cell(self):
        return self._cell
    @cell.setter
    def cell(self, cell):
        from ase.cell import Cell
        cell = Cell.ascell(cell)
        self._cell[:] = cell
    def set_cell(self, cell, scale_atoms=False):
        """Set unit cell vectors.

        Parameters:

        Examples:

        """
        from ase.cell import Cell
        from ase.geometry.cell import complete_cell

        cell = Cell.new(cell)
        oldcell = Cell(self.cell)
        self.cell = cell
        if scale_atoms:
            M = np.linalg.solve(oldcell.complete(), cell.complete())
            for ba in self.batoms.values():
                ba.set_positions(np.dot(ba.get_positions(), M))
    @property
    def pbc(self):
        return self.get_pbc()
    @pbc.setter
    def pbc(self, pbc):
        self.set_pbc(pbc)
    def get_pbc(self):
        return list(self.coll.blase.pbc)
    def set_pbc(self, pbc):
        if isinstance(pbc, bool):
            pbc = [pbc]*3
        self.coll.blase.pbc = pbc
    @property
    def boundary(self):
        return self.get_boundary()
    @boundary.setter
    def boundary(self, boundary):
        self.update_boundary(boundary)
    def get_boundary(self):
        boundary = np.array(self.coll.blase.boundary)
        return boundary.reshape(3, -1)
    def update_boundary(self, boundary = None):
        """
        >>> from blase.batoms import Batoms
        >>> from ase.io import read
        >>> atoms = read('docs/source/_static/datas/tio2.cif')
        >>> tio2 = Batoms(label = 'tio2', atoms = atoms, model_type = '2', polyhedra_dict = {'Ti': ['O']}, color_style="VESTA")
        >>> tio2.boundary = 0.5
        """
        from ase.cell import Cell
        tstart = time()
        self.clean_blase_objects('boundary')
        if boundary is not None:
            if isinstance(boundary, (int, float)):
                boundary = np.array([-boundary, 1 + boundary]*3)
            elif len(boundary) == 3:
                if isinstance(boundary[0], (int, float)):
                    boundary = np.array([[-boundary[0], 1 + boundary[0]],
                                      [-boundary[1], 1 + boundary[1]],
                                      [-boundary[2], 1 + boundary[2]]])
                elif len(boundary[0]) == 2:
                    boundary = np.array(boundary)
            else:
                raise Exception('Wrong boundary setting!')
            self.coll.blase.boundary = boundary[:].flatten()
        boundary = self.boundary
        atoms_skin = Atoms()
        atoms_skin.info['species'] = []
        if self.atoms.pbc.any():
            for species, batom in self.batoms.items():
                positions1, positions2 = search_boundary(batom.local_positions, self.cell, boundary)
                positions1 = positions1 + batom.location
                ba = Batom(self.label, '%s_boundary'%(species), positions1, scale = batom.scale, 
                            segments = self.segments, shape = self.shape, material=batom.material)
                self.coll.children['%s_boundary'%self.label].objects.link(ba.batom)
                atoms_skin = atoms_skin + Atoms('%s'%batom.element*len(positions2), positions2)
                atoms_skin.info['species'].extend([species]*len(positions2))
            # print('search boundary: {0:10.2f} s'.format(time() - tstart))
            tstart = time()
            # find skin
            atoms, n1, n2, n3 = self.get_atoms_with_boundary()
            skin = list(range(len(atoms), len(atoms) + len(atoms_skin)))
            atoms = atoms + atoms_skin
            atoms.info['species'].extend(atoms_skin.info['species'])
            self.bondlist = build_bondlists(atoms, self.bondsetting.data)
            atoms_skin = search_skin(atoms, self.bondsetting.data, self.bondlist, skin)
            self.update_skin(atoms_skin)
            # print('search skin: {0:10.2f} s'.format(time() - tstart))
    def update_skin(self, atoms):
        # print(atoms)
        self.clean_blase_objects('skin')
        if len(atoms) == 0: return
        tstart = time()
        specieslist = set(atoms.info['species'])
        for species in specieslist:
            ind = [i for i, x in enumerate(atoms.info['species']) if x == species]
            positions = atoms.positions[ind]
            ba = Batom(self.label, '%s_skin'%(species), positions, scale = self.batoms[species].scale,
                        segments=self.segments, shape=self.shape)
            self.coll.children['%s_skin'%self.label].objects.link(ba.batom)
        # print('update skin: {0:10.2f} s'.format(time() - tstart))

    @property
    def model_type(self):
        return self.get_model_type()
    @model_type.setter
    def model_type(self, model_type):
        self.set_model_type(model_type)
    def get_model_type(self):
        return np.array(self.coll.blase.model_type)
    def set_model_type(self, model_type):
        self.coll.blase.model_type = str(model_type)
        self.draw()
    @property
    def show_unit_cell(self):
        return self.get_show_unit_cell()
    @show_unit_cell.setter
    def show_unit_cell(self, show_unit_cell):
        self.set_show_unit_cell(show_unit_cell)
    def get_show_unit_cell(self):
        return self.coll.blase.show_unit_cell
    def set_show_unit_cell(self, show_unit_cell):
        self.coll.blase.show_unit_cell = show_unit_cell
        self.draw_cell()
    @property
    def atoms(self):
        return self.get_atoms()
    def get_atoms(self):
        return self.batoms2atoms(self.batoms)
    @property
    def atoms_boundary(self):
        return self.get_atoms_boundary()
    def get_atoms_boundary(self):
        return self.batoms2atoms(self.batoms_boundary)
    @property
    def atoms_skin(self):
        return self.get_atoms_skin()
    def get_atoms_skin(self):
        return self.batoms2atoms(self.batoms_skin)
    def get_atoms_with_boundary(self):
        """
        build ASE atoms from batoms dict.
        """
        atoms = self.atoms
        n1 = len(atoms)
        atoms_boundary = self.atoms_boundary
        n2 = len(atoms_boundary)
        atoms_skin = self.atoms_skin
        n3 = len(atoms_skin)
        species = atoms.info['species'] + atoms_boundary.info['species'] + atoms_skin.info['species']
        atoms = atoms + atoms_boundary + atoms_skin
        atoms.info['species'] = species
        atoms.pbc = False
        return atoms, n1, n2, n3
    @property
    def positions(self):
        return self.get_positions()
    def get_positions(self):
        return self.atoms.positions
    def get_scaled_positions(self):
        return self.atoms.get_scaled_positions()
    @property
    def species(self):
        return self.get_species()
    def get_species(self):
        """
        build species from collection.
        """
        species = []
        for ba in self.coll_atom.objects:
            species.append(ba.species)
        return species
    @property
    def batoms(self):
        return self.get_batoms()
    def get_batoms(self):
        batoms = {}
        for ba in self.coll_atom.objects:
            batoms[ba.species] = Batom(ba.name)
        return batoms
    @property
    def batoms_boundary(self):
        return self.get_batoms_boundary()
    def get_batoms_boundary(self):
        batoms_boundary = {}
        for ba in self.coll_boundary.objects:
            batoms_boundary[ba.species] = Batom(ba.name)
        return batoms_boundary
    @property
    def batoms_skin(self):
        return self.get_batoms_skin()
    def get_batoms_skin(self):
        batoms_skin = {}
        for ba in self.coll_skin.objects:
            batoms_skin[ba.species] = Batom(ba.name)
        return batoms_skin
    def batoms2atoms(self, batoms, local = False):
        object_mode()
        atoms = Atoms()
        species_list = []
        symbols = []
        positions = []
        for species, batom in batoms.items():
            if species[-9:] == '_boundary': species = species[0:-9]
            if species[-5:] == '_skin': species = species[0:-5]
            species_list.extend([species]*len(batom))
            symbol = [batom.element]*len(batom)
            symbols.extend(symbol)
            if local:
                positions.extend(batom.local_positions)
            else:
                positions.extend(batom.positions)
        atoms = Atoms(symbols, positions, cell = self.cell, pbc = self.pbc)
        atoms.info['species'] = species_list
        return atoms
    def draw_constraints(self):
        """
        """
        #
        constr = self.atoms.constraints
        self.constrainatoms = []
        for c in constr:
            if isinstance(c, FixAtoms):
                for n, i in enumerate(c.index):
                    self.constrainatoms += [i]
    
    def highlight_atoms(self, indexs, shape = 'sphere', radius_scale=1.4,
                           color=(0.0, 1.0, 0.0), transmit=0.6):
        """
        """
        object_mode()
        for index in indexs:
            loc = self.positions[index]
            ele = self.symbols[index]
            radii = radius_scale * self.atom_kinds[ele]['radius']
            if shape == 'cube':
                bpy.ops.mesh.primitive_cube_add(location=loc, size=radii*2)
            else:
                bpy.ops.mesh.primitive_uv_sphere_add(location=loc, radius=radii)
            ball = bpy.context.view_layer.objects.active
            bpy.ops.object.shade_smooth()
            ball.data.materials.append(material)
            ball.show_transparent = True
            coll_highlight.objects.link(ball)
    def load_frames(self, images = None):
        """

        images: list
            list of atoms. All atoms show have same species and length.
            
        >>> from ase.io import read
        >>> from blase import Batoms
        >>> images = read('h2o-animation.xyz', index = ':')
        >>> h2o = Batoms(label = 'h2o', atoms = images)
        >>> h2o.load_frames()
        >>> h2o.render(animation = True)

        """
        if not images:
            images = self.images
        if len(self.atoms) != len(images[0]):
            raise Exception("Number of atoms %s is not equal to %s."%(len(self.atoms), len(images[0])))
        atoms = images[0]
        if 'species' not in atoms.info:
            atoms.info['species'] = atoms.get_chemical_symbols()
        positions = np.array([atoms.positions for atoms in images])
        for species, ba in self.batoms.items():
            index = [atom.index for atom in atoms if atoms.info['species'][atom.index] == species]
            ba.load_frames(positions[:, index])
    
    def calc_camera_data(self, canvas, canvas1, direction = (0, 0, 1)):
        """
        """
        from scipy.spatial.transform import Rotation as R
        camera_target = np.mean(canvas, axis=0)
        camera_data = {}
        width = canvas1[1, 0] - canvas1[0, 0]
        height = canvas1[1, 1] - canvas1[0, 1]
        ortho_scale = max(width, height)
        #
        direction = direction/np.linalg.norm(direction)
        location = camera_target + direction*20
        camera_data = {'camera_loc': location, 'camera_target': camera_target,
                        'ortho_scale': ortho_scale, 'ratio': height/width}
        return camera_data
    def calc_bond_data(self, atoms, bondlists):
        """
        """
        from blase.tools import get_bond_kind
        from ase.data import chemical_symbols, covalent_radii
        batoms = self.batoms
        positions = atoms.positions
        chemical_symbols = np.array(chemical_symbols)
        if 'species' not in atoms.info:
            atoms.info['species'] = atoms.get_chemical_symbols()
        speciesarray = np.array(atoms.info['species'])
        elementarray = np.array(atoms.get_chemical_symbols())
        radiusarray = np.array([covalent_radii[chemical_symbols == em][0] for em in elementarray])
        scalearray = np.array([batoms[sp].scale for sp in speciesarray])
        specieslist = list(set(atoms.info['species']))
        bond_kinds = {}
        if len(bondlists) == 0:
            self.bond_kinds = bond_kinds
            return
        for spi in specieslist:
            bondlists1 = bondlists[speciesarray[bondlists[:, 0]] == spi]
            if len(bondlists1) == 0: continue
            emi = spi.split('_')[0]
            emj = elementarray[bondlists1[:, 1]]
            bond_kind = get_bond_kind(emi, color_style=self.color_style)
            if spi not in bond_kinds:
                bond_kinds[spi] = bond_kind
            offset = bondlists1[:, 2:5]
            R = np.dot(offset, atoms.cell)
            vec = positions[bondlists1[:, 0]] - (positions[bondlists1[:, 1]] + R)
            length = np.linalg.norm(vec, axis = 1)
            nvec = vec/length[:, None]
            radius1 = covalent_radii[chemical_symbols ==emi][0]
            radius2 = radiusarray[bondlists1[:, 1]]
            scale2 = scalearray[bondlists1[:, 1]]
            pos = [positions[bondlists1[:, 0]] - nvec*radius1*batoms[spi].scale*0.5,
                    positions[bondlists1[:, 1]] + R + (nvec.T*radius2).T*scale2*0.5]
            center0 = (pos[0] + pos[1])/2.0
            vec = pos[0] - pos[1]
            length = np.linalg.norm(vec, axis = 1)
            nvec = vec/length[:, None]
            nvec = nvec + 1e-8
            # verts, faces
            v1 = nvec + np.array([1.2323, 0.493749, 0.5604937284])
            tempv = np.einsum("ij, ij->i", v1, nvec)
            v11 = v1 - (nvec.T*tempv).T
            templengh = np.linalg.norm(v11, axis = 1)
            v11 = v11/templengh[:, None]/2.828427
            tempv = np.cross(nvec, v11)
            v22 = (tempv.T*(length*length)).T
            #
            center = (center0 + pos[0])/2.0
            bond_kinds[spi]['centers'] = center
            bond_kinds[spi]['lengths'] = length/4.0
            bond_kinds[spi]['normals'] = nvec
            bond_kinds[spi]['verts'] = center + v11
            bond_kinds[spi]['verts'] = np.append(bond_kinds[spi]['verts'], center - v11, axis = 0)
            bond_kinds[spi]['verts'] = np.append(bond_kinds[spi]['verts'], center + v22, axis = 0)
            bond_kinds[spi]['verts'] = np.append(bond_kinds[spi]['verts'], center - v22, axis = 0)
        self.bond_kinds = bond_kinds
    

    def get_distances(self, species1, i, species2, indices, mic=False):
        """
        Return distances of atom No.i with a list of atoms.

        Use mic=True to use the Minimum Image Convention.

        >>> h2o.get_distances('O', 0, 'H', [0, 1])
        """
        from ase.geometry import get_distances

        p1 = self.batoms[species1][i]
        p2 = self.batoms[species2][indices]
        cell = None
        pbc = None
        if mic:
            cell = self.cell
            pbc = self.pbc
        D, D_len = get_distances(p1, p2, cell=cell, pbc=pbc)
        D_len.shape = (-1,)
        return D_len
    def get_angle(self, species1, i1, species2, i2, species3, i3, mic=False):
        """
        Get angle in degrees between the vectors i2->i1 and
        i2->i3.
        Use mic=True to use the Minimum Image Convention and calculate the
        angle across periodic boundaries.

        >>> h2o.get_angle('H', 0, 'O', 0, 'H', 1)

        """
        from ase.geometry import get_angles
        p1 = self.batoms[species1][i1]
        p2 = self.batoms[species2][i2]
        p3 = self.batoms[species3][i3]
        v12 = p1 - p2
        v32 = p3 - p2
        cell = None
        pbc = None
        if mic:
            cell = self.cell
            pbc = self.pbc
        return get_angles([v12], [v32], cell=cell, pbc=pbc)
