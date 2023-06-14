# Welcome to PyBEAM!
PyBEAM (Bayesian Evidence Accumulation Models) is a Python package designed to rapidly fit two-threshold, binary choice models to choice-RT data using Bayesian inference methods. For a full description of its design, see the publication (https://psyarxiv.com/ax36b/). For access to the package code and other files, see the PyBEAM github (https://github.com/murrowma/pybeam/). To learn how to install and use PyBEAM, see the PyBEAM documentation (https://pybeam-documentation.readthedocs.io/en/latest/index.html) and tutorials/example files located here.

# Installation

## System requirements and package dependencies
- Python 3.9 or higher
- PyMC (v5) and its dependencies (automatically installed when PyMC is installed)
- numpy, matplotlib, and ArviZ (installed automatically with PyMC)
- Cython

## Installing PyBEAM
To install PyBEAM, you must first install PyMC. To do so, follow the instructions located on PyMC’s website. The url for this as of April 17, 2023 is:

  https://www.pymc.io/projects/docs/en/latest/installation.html

Once you have installed PyMC, you must next install the Python package Cython and a C compilier. PyBEAM utilizes Cython to improve program speed, requiring files to be compilied. Cython is installed via pip or conda like other Python packages, using either:
  
  pip install Cython
  conda install Cython

More detailed instructions for installing Cython are available at the following link if needed:

  https://cython.readthedocs.io/en/stable/src/quickstart/install.html

After Cython is installed, download a compilier. For Mac, this is part of Xtools, which can be downloaded from terminal or the app store. If it is not already on your Mac, you will be prompted automatically to install it when you install PyBEAM. To download directly from terminal, open a terminal window and input the following command:

  xcode-select --install

For Windows, a compilier must be downloaded separately. A couple options are available to do this. The first is to download Window’s Visual Studio program. Instructions for this are at the following url (as of December 15, 2022):

  https://github.com/cython/cython/wiki/CythonExtensionsOnWindows

The second option is to install MinGW (NOT MinGW-64). Though this is not recommended, instructions for this are located on Cython’s docs, located at this URL (as of December 15, 2022):

  https://cython.readthedocs.io/en/latest/src/tutorial/appendix.html?highlight=compile%20windows

We can now install PyBEAM. First, download and unzip pybeam.zip. Then, open one of conosle, command prompt, or anaconda prompt (depending on your system), and navigate to the pybeam directory. In this directory there should be two files: the folder pybeam, and a setup.py file. Once there, run the following line of code (INCLUDING the period):

  pip install -e .

Once you have done this, PyBEAM is ready to be used!

## Learning to use PyBEAM
PyBEAM contains two submodules: precoded and custom. The precoded submodule provides many precoded models from the literature which should be sufficient for most needs. The second submodule, custom, provides tools for users to create their own models outside the scope of the precoded set. Tutorials for using the models are provided in Jupyter notebook form under the Precoded tutorials and Custom tutorials tabs. For full documentation of PyBEAM’s functions, see tabs Precoded functions and Custom functions. Futher examples are provided in folder “Examples” on the PyBEAM github.

## Custom models
Folder “custom_model_template.zip” on the PyBEAM github contains the custom model template files. See the Custom tutorials tab for instructions on how to create your own model.
