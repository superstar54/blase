import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="blendase",
    version="1.0.0",
    description="Drawing and rendering atoms, molecule and crystal using Blender.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/superstar54/blendase",
    author="Xing Wang",
    author_email="xingwang1991@gmail.com",
    license="GPL",
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        "Programming Language :: Python :: 3.9",
    ],
    packages=["blendase"],
    #include_package_data=True,
    package_data = {
    'blendase': ['docs/source/_static/*png'],
    },
    install_requires=["ase", "numpy", "scipy", "scikit-image"],
    python_requires='>=3.9',
)
