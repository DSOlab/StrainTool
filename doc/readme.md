# General

# Installation

## Prerequisites

The package (actually called `pystrain`) has been sucesefuly installed under
both python2.7 and python3.6.
Apart from a/the python interpreter, the package will need Numpy and Scipy.

The main program tha computes the strain tensors is located at `StrainTool/bin`,
named `StrainTensor.py`. Before you go ahead and use it though, you need to
install the `pystrain` package.

## `pystrain` package installation
Go to the `StrainTool/pystrain` directory and run the command
`python setup.py install`. You will probably need to do this with root 
privileges.

## `StrainTensor.py`

This is the program that actualy does the computations, found under `StrainTool/bin`.
It may be needed to make it executable before you run it (in a UNIX-like OS, 
the following should do `chmod +x StrainTensor.py`).
