.. _download_and_install:

============
Installation
============

Requirements
============
* Blender_ 2.93 or newer
* Python_ 3.9 or newer
* ASE_ 

Optional:

* scikit-image_ for isosurface
* pymatgen_ for Pymatgen package

.. _Blender: https://www.blender.org/
.. _Python: https://www.python.org/
.. _ASE: https://wiki.fysik.dtu.dk/ase/index.html
.. _Pymatgen: https://pymatgen.org/
.. _scikit-image: https://scikit-image.org/


.. index:: pip
.. _pip installation:




Install ASE
======================

On Windows, suggest to install Python with ``Anaconda``, https://docs.anaconda.com/anaconda/install/windows/

First install ASE on your computer. On Windows, open Anaconda Prompt, on Linux, open a terminal, and run::
    
    pip3 install --upgrade ase

Then install ASE inside Blender. On Linux, go to your Blender python directory, e.g. ``blender-2.93-linux-x64/2.93/python/bin``, install pip_::
    
    $ ./python3.9 -m ensurepip
    
Install ASE_ and scikit-image_ inside Blender::

    $ ./pip3 install --upgrade ase
    
    $ ./pip3 install scikit-image

On Windows, start Blender, open a Python console and run the following code::

    import sys
    import os
    import subprocess
    path = os.path.join(sys.prefix, 'bin', 'python.exe')
    subprocess.call([path, "-m", "ensurepip"])
    subprocess.call([path, "-m", "pip", "install", "--upgrade", "pip"])
    subprocess.call([path, "-m", "pip", "install", "ase"])
    subprocess.call([path, "-m", "pip", "install", "scikit-image"])
 
.. note::

   You could avoid install inside Blender by setting bl to use system python package::

    export BLENDER_COMMAND='blender --python-use-system-env'


Install Blase
========================
You can get the source from https://github.com/superstar54/blase.

:Zip-file:

    You download the source as a `zip-file`.

:Git clone:

    Alternatively, you can get the source for the latest stable release like this::

        $ git clone --depth 1 https://github.com/superstar54/blase.git

- Extract the file, rename the file to ``blase``, and move it to ``blender-2.93.4-linux-x64/2.93/scripts/addons/``. 

- Enable the addon in the Preferences setting. Please vist here to learn how to enable an addon. 
https://docs.blender.org/manual/en/latest/editors/preferences/addons.html. Or, you can open a Blender Python console, and run the following code to enable the blase addon::

    import addon_utils
    addon_utils.enable('blase', default_set=True)


Environment variables
=====================

First, we need to get the full path to the blase package. Start Blender, in the python console, run::

    import blase
    blase.__path__

On Linux, set these permanently in your :file:`~/.bashrc` file::

    export BLENDER_COMMAND='blender --python-use-system-env'
    export PYTHONPATH=<path-to-blase-package>:$PYTHONPATH
    export PATH=<path-to-blase-package>bin/:$PATH
    export BLASE_PATH=<path-to-blase-package>/bin/

.. note::

   On windows, you can edit the system environment variables.



Test your installation
======================

Start Blender, in the python console, run:

>>> from blase.batoms import Batoms
>>> h2o = Batoms({'O': [[0, 0, 0.40]], 'H': [[0, -0.76, -0.2], [0, 0.76, -0.2]]})


.. image:: _static/batoms-h2o.png
   :width: 3cm
   
If you saw a water molecule, congratulations!


Install Pymatgen
================================

Rename you blender python folder (``blender-2.93-linux-x64/2.93/python``) to ``_python``. Create a virtual environment for your blender using conda::

    conda create --prefix $Path_to_blener/blender-2.93.4-linux-x64/2.93/python python=3.9.2


On Linux, go to the new python directory, e.g. ``blender-2.93-linux-x64/2.93/python/bin``, and install ASE_,  scikit-image_ and Pymatgen_ inside Blender::

    $ ./pip3 install --upgrade ase
    
    $ ./pip3 install scikit-image

    $ ./pip3 install pymatgen
