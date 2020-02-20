from setuptools import setup

setup(name = 'blase',
      version='0.1',
      description='A blender interface to ASE',
      url='http://github.com/superstar54/blase',
      author='Xing Wang',
      author_email='xingwang1991@gmail.com',
      license='MIT',
      platforms=['linux'],
      packages=['blase'],
      long_description='''Python module for drawing and rendering ase atoms objects using blender.''',
      install_requires=[
          "blender"
      	  "ase"
          "numpy",
          "scipy",
          "skimage"
],)