import numpy as np
import pickle
import os
import json
from ase.data import covalent_radii, atomic_numbers, chemical_symbols
from ase.data.colors import jmol_colors
import pprint



def get_bondpairs(atoms, cutoff=1.0, rmbonds = []):
    """
    Get all pairs of bonding atoms
    rmbonds
    """
    from ase.data import covalent_radii
    from ase.neighborlist import NeighborList
    cutoffs = cutoff * covalent_radii[atoms.numbers]
    nl = NeighborList(cutoffs=cutoffs, self_interaction=False)
    nl.update(atoms)
    bondpairs = []
    natoms = len(atoms)
    for a in range(natoms):
        indices, offsets = nl.get_neighbors(a)
        # print(indices, offsets)
        for a2, offset in zip(indices, offsets):
            flag = True
            for rmpair in rmbonds:
                if atoms[a].symbol == rmpair[0] and atoms[a2].symbol == rmpair[1] \
                  or atoms[a].symbol == rmpair[1] and atoms[a2].symbol == rmpair[0]:
                    flag = False
            # print(a, a2, flag)
            if flag:
                bondpairs.extend([([a, a2], offset)])
    return bondpairs

def write_blender(atoms, display = False, queue = None, **kwargs):
    with open('blase.inp', 'wb') as f:
        pickle.dump([atoms, kwargs], f)
    #
    blender_cmd = os.environ['BLENDER_COMMAND']
    blase_path = os.environ['BLASE_PATH']
    blase_cmd = blase_path + '/bin/run-blase.py'
    if display:
        cmd = blender_cmd + ' -P ' + blase_cmd
    elif queue == 'SLURM':
        cmd = 'srun -n $SLURM_NTASKS ' +  blender_cmd + ' -b ' + ' -P ' + blase_cmd
    else:
        cmd = blender_cmd + ' -b ' + ' -P ' + blase_cmd
    # print(cmd)
    errcode = os.system(cmd)
    # if errcode != 0:
    #     raise OSError('Command ' + cmd +
    #                   ' failed with error code %d' % errcode)

# def get_atom_kinds(atoms, props):
    # return kinds
def get_atom_kinds(atoms, props = {}):
    # symbols = atoms.symbols
    # formula = atoms.symbols.formula
    # atom_kinds = formula.count()
    if hasattr(atoms, 'kinds'):
        kinds = list(set(atoms.kinds))
    else:
        atoms.kinds = atoms.get_chemical_symbols()
        kinds = list(set(atoms.kinds))
    print(kinds)
    atom_kinds = {}
    for kind in kinds:
        atom_kinds[kind] = {}
        element = kind.split('_')[0]
        number = chemical_symbols.index(element)
        inds = [atom.index for atom in atoms if atoms.kinds[atom.index]==kind]
        color = jmol_colors[number]
        radius = covalent_radii[number]
        atom_kinds[kind]['element'] = element
        atom_kinds[kind]['positions'] = atoms[inds].positions
        atom_kinds[kind]['number'] = number
        atom_kinds[kind]['color'] = color
        atom_kinds[kind]['transmit'] = 1.0
        atom_kinds[kind]['radius'] = radius
        atom_kinds[kind]['balltype'] = None
        if props:
            if kind in props.keys():
                for prop, value in props[kind].items():
                    atom_kinds[kind][prop] = value
    return atom_kinds
def get_bond_kinds(atoms, bondlist):
    '''
    Build faces for instancing bonds.
    The radius of bonds is determined by nbins.
    mesh.from_pydata(vertices, [], faces)
    '''
    bond_kinds = {}
    for bond in bondlist:
        inds, offset = bond
        R = np.dot(offset, atoms.cell)
        center0 = atoms.positions[inds[0]] + atoms.positions[inds[1]]
        # kinds = [atoms[ind].symbol for ind in [a, b]]
        i = 0
        for ind in inds:
            kind = atoms[ind].symbol
            if kind not in bond_kinds.keys():
                lengths = []
                centers = []
                normals = []
                bond_kinds[kind] = {'lengths': lengths, 'centers': centers, 'normals': normals}
                number = chemical_symbols.index(kind)
                color = jmol_colors[number]
                radius = covalent_radii[number]
                bond_kinds[kind]['number'] = number
                bond_kinds[kind]['color'] = color
                bond_kinds[kind]['transmit'] = 1.0
            #
            center = (0.5 * (center0 + (-1)**i*R) + atoms[ind].position)/2.0
            if atoms.positions[ind][2] > center[2]:
                vec = atoms.positions[ind] - center
            else:
                vec = center - atoms.positions[ind]
            length = np.linalg.norm(vec)
            nvec = vec/length
            bond_kinds[kind]['lengths'].append(length)
            bond_kinds[kind]['centers'].append(center)
            bond_kinds[kind]['normals'].append(nvec)
            i += 1
    # pprint.pprint(bond_kinds)
    return bond_kinds

def euler_from_vector(normal, s = 'zxy'):
    from scipy.spatial.transform import Rotation as R
    normal = normal/np.linalg.norm(normal)
    vec = np.cross([0.0000014159, 0.000001951, 1], normal)
    vec = vec/np.linalg.norm(vec)
    # print(vec)
    # ang = np.arcsin(np.linalg.norm(vec))
    ang = np.arccos(normal[2])
    vec = -1*ang*vec
    # print(vec)
    r = R.from_rotvec(vec)
    euler = r.as_euler(s)
    return euler

ang = euler_from_vector(np.array([1, 1, 2]))
print(ang)



def getEquidistantPoints(p1, p2, n):
    return zip(np.linspace(p1[0], p2[0], n+1), np.linspace(p1[1], p2[1], n+1), np.linspace(p1[2], p2[2], n+1))

