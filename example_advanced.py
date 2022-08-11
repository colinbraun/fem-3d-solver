# Import matplotlib
import matplotlib.pyplot as plt
# Import numpy
import numpy as np
# Import Waveguide3D
from waveport.waveport import Waveguide3D

# Construct the Waveguide3D object. This loads the mesh into the object and initializes variables
# This particular geometry is oriented along the z-axis (this was chosen when it was created in Cubit)
# The volume names, volume permittivities, and PEC name
vn, vp, pn = ["TetrahedronsVacuum", "TetrahedronsSubstrate"], [1, 4.5], "PECWalls"
# The Absorbing Boundary Condition walls name (if it exists)
abcn = "ABC"
# The input port surface names, the input port boundary name, and the input port permittivities
# There are two input port surfaces for the microstrip line. The substrate part and the vacuum/air part.
ipn, ipbn, ipp = ["InputPortVacuum", "InputPortSubstrate"], "InPortPEC", [1, 4.5]
# The output port surface names, the output port boundary name, and the output port permittivities
# There are two output port surfaces for the microstrip line. The substrate part and the vacuum/air part.
opn, opbn, opp = ["OutputPortVacuum", "OutputPortSubstrate"], "OutPortPEC", [1, 4.5]
# Integration lines. These are optional, but if included can guarantee the polarity of the fields in a simulation.
# Not doing this will give correct field results, but the input and output port field profiles used in the
# S-parameter calculations may be opposite of each other, causing the phase of S21 to be incorrect.
# This idea is identical to the integration lines you have to choose in HFSS.
p1ip, p2ip = np.array([0, 0.0008]), np.array([0, 0])
p1op, p2op = np.array([0, 0.0008]), np.array([0, 0])
# The chosen operating frequency. k0 = 2*pi*f/c to convert to/from frequency.
k0 = 4
# Create the Waveguide object with those specifications
waveguide = Waveguide3D("meshes/microstrip_line_44000tets_pec_walls_example.inp", k0, vn, vp, pn, abcn, ipn, ipbn, ipp, opn, opbn, opp, p1ip, p2ip, p1op, p2op)

# Everything else is now the same as the simple example.
# Solve the FEM problem. This constructs the matrix equation and solves it by taking a matrix inverse
waveguide.solve()

# Determine how long the geometry is along the z-axis. z_max and z_min are attributes of the Waveguide3D object
z_length = waveguide.z_max - waveguide.z_min

# Produce some field results in the XY plane offset halfway along the geometry. This produces a matplotlib figure.
# On Bell, the interactive mode of matplotlib does not seem to work (it won't pop up an figure to interact with)
waveguide.plot_fields(plane="xy", offset=z_length/2, phase=0)

# matplotlib still works without interactive mode, but you will need to save the the results to view them
plt.savefig("fields_plot_advanced.png")
