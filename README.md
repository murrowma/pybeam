# Welcome to PyBEAM!
PyBEAM (Bayesian Evidence Accumulation Models) is a Python package designed to rapidly fit two-boundary, binary choice models to choice-RT data using Bayesian inference methods. For a full description of its design, see the publication (future link here). To learn how to use PyBEAM, see the tutorials' directory for a step by step introduction.

# Installation

To install PyBEAM, you must first install a C compilier. PyBEAM utilizes Cython to improve program speed, requiring a files to be compilied. For Mac, this is part of the xcode library. If it is not already on your Mac, you will be prompted automatically to install it when you install PyBEAM. For Windows, a compilier must be downloaded separately. A couple options are available to do this. The first is to download Window's Visual Studio program. Instructions for this are at the following url (as of January 30, 2022):

    https://github.com/cython/cython/wiki/CythonExtensionsOnWindows
    
The second option is to install MinGW (NOT MinGW-64). Instructions for this are located on Cython's docs, located at this URL (as of January 30, 2022):

    https://cython.readthedocs.io/en/latest/src/tutorial/appendix.html?highlight=compile%20windows
    
After C compiliers have been installed, you must next install Cython. This can be pip installed like other Python packages. Instructions for this are at the follwoing link:

    https://cython.readthedocs.io/en/stable/src/quickstart/install.html
    
Once you have downloaded and linked a compilier to Python, you must now download PyMC3. Instructions for this are located on PyMC3's website at the url below (as of January 30, 2022):

    https://docs.pymc.io/en/v3/
    
Once you have installed PyMC3, download and unzip pybeam.zip. Then, open one of conosle, command prompt, or anacond prompt (depending on your system), and navigate to the pybeam directory. In this directory there should be two files: the folder pybeam, and a setup.py file. Once there, run the following line of code (INCLUDING the period):

    pip install .
    
Once you have done this, PyBEAM is ready to be used!

# Learning to use PyBEAM

PyBEAM contains two sub-modules: default and custom. The default sub-module provides many pre-coded models from the literature which should be sufficient for most needs. The second sub-module, custom, provides tools for users to create their own models outside the scope of the pre-coded set.

Jupyter notebook tutorials which introduce both sub-modules are provided here. The folder "default_tutorials" provides notebooks for the default sub-module, while "custom_tutorials" provides tutorials for the cutsom sub-module.

# Custom models

Folder "custom_model_template.zip" contains the custom model template files. See "Custom_Tutorial1_creating_a_custom_model" for directions on how to use these files to create your own custom model.
