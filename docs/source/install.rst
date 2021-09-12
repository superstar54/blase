.. _download_and_install:

============
Installation
============

Requirements
============
* Blender_ 2.90 or newer
* Python_ 3.7 or newer
* ASE_ 

Optional:

* scikit-image_ for isosurface

.. _Blender: https://www.blender.org/
.. _Python: https://www.python.org/
.. _ASE: https://wiki.fysik.dtu.dk/ase/index.html
.. _scikit-image: https://scikit-image.org/


.. index:: pip
.. _pip installation:



Install Python
=====================

On Windows, suggest to install Python with ``Anaconda``, https://docs.anaconda.com/anaconda/install/windows/


Install ASE
======================

.. highlight:: bash

First install ASE on your computer. On Windows, open Anaconda Prompt, on Linux, open a terminal, and run::
    
    pip3 install --upgrade ase

Then install ASE inside Blender. On Linux, go to your Blender python directory, e.g. ``blender-2.92/2.92/python/bin``, install pip_::
    
    $ ./python3 -m ensurepip
    
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

   You could avoild install inside Blender by setting bl to use system python package::

    export BLENDER_COMMAND='blender --python-use-system-env'


Install Blase
========================
You can get the source from https://github.com/superstar54/blase.

:Zip-file:

    You download the source as a `zip-file`.

:Git clone:

    Alternatively, you can get the source for the latest stable release like this::

        $ git clone --depth 1 https://github.com/superstar54/blase.git
    And then compress the blase folder to a `zip-file`.

Then install blase as a addon in Blender. Please vist here to learn how to install a addon with zip files. 
https://docs.blender.org/manual/en/latest/editors/preferences/addons.html. 


Don't forget to enable the addon. You can enable in the Preferences setting. Or, open a Blender Python console, and run::

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

Before running the tests, make sure you have set your :envvar:`PATH`
environment variable correctly as described in the relevant section above.
Run the tests like this::

    $ blase  

If you see a water molecule in a Blender window, then your installation is successful.

