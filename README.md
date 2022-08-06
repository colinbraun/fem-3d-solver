An initial README file. Will be updated at some point.

# 3D FEM Solver

## Overview
This software provides a semi-general-purpose FEM simulator for 3D geometries. It takes a .inp file (usually generated in Coreform Cubit) and performs a simulation on it. The code is currently set up to have 1 input port and 1 output port.

# Explanation of How it works
There are several steps that must be done in order to perform a simulation.
1. The mesh data from the .inp file must be loaded and manipulated such that the salient data is available to the rest of the program.
2. A matrix equation must be constructed.
3. The matrix equation must be solved.
4. Meaningful results must be generated.

Each of these sections will be addressed next.

## Mesh Loading and Data Structuring
To understand this section best, you should probably have a .inp file to look at as an example, as well as Coreform Cubit pulled up. It is assumed we are working with tetrahedrons.

For a 3D mesh of tetrahedrons created in Cubit and exported to a .inp file, Cubit will create a section in the file containing the locations of each node, assigning each a unique number to identify them (essentially a global node number). This section is labeled "ALLNODES" in the .inp file (not sure if you can change this or not). This section will looks something like this in the .ipn file:

********************************** N O D E S **********************************
*NODE, NSET=ALLNODES
       1,    5.000000e-01,   -2.500000e-01,    1.750000e+00
       2,    5.000000e-01,    2.500000e-01,    1.750000e+00
       3,   -5.000000e-01,    2.500000e-01,    1.750000e+00
       4,   -5.000000e-01,   -2.500000e-01,    1.750000e+00
       5,    5.000000e-01,    2.500000e-01,    0.000000e+00

One can have Cubit provide more information about the mesh than this by creating what Cubit calls "Blocks". You can see this in the GUI. These Blocks contain information about which node numbers (call them global node numbers if you prefer) make up a particular surface, volume, curve/edge, etc. This information is given in the .inp always in terms of global node numbers.

Suppose we are working with a cuboid as our volume (perhaps it is a rectangular waveguide). We would create a Block for the cuboid volume in Cubit. When we export the .inp file, we will see the following kind of information:

********************************** E L E M E N T S ****************************
*ELEMENT, TYPE=C3D4, ELSET=Tetrahedrons
       1,     824,     825,     826,     827
       2,     827,     829,     828,     830
       3,     827,     829,     830,     831
       4,     825,     827,     830,     831
       5,     832,     833,     834,     835

These are just the first 5 entries of the Block. Each entry describes the 4 nodes of a tetrahedron that makes up part of the voluem. Notice they contain the global node number from the ALLNODES set earlier. It is also worth nothing that Blocks can be named. This one is named "Tetrahedrons", as visible above. This was a choice, and not automatic. The default names for blocks are usually things like "EB1", "EB2", and so on. It is just a convenience to name them something meaningful.

Similarly, if we created a Block containing information about a surface, each entry would contain 3 global node numbers that make up a triangle on the surface.

### Current Mesh Loading and Data Structuring Implementation
All of the mesh loading and data structuring is done in util.py. A function called `load_mesh_block()` will take a string containing the name of the Block and return a 2D array containing the data.


