from batoms.batoms import Batoms
import pickle

def run():
    with open('.batoms.inp', 'rb') as f:
        atoms, batoms_input, render_input = pickle.load(f)
        if 'label' not in batoms_input:
            batoms_input['label'] = 'batoms'
        batoms = Batoms(atoms = atoms, **batoms_input)
        batoms.render.run(**render_input)

if __name__ == "__main__":
    run()
