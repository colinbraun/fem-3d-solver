# Import matplotlib
import matplotlib.pyplot as plt
# Import Waveguide3D
from waveport.waveport import Waveguide3D

# Construct the Waveguide3D object. This loads the mesh into the object and initializes variables
# This particular geometry is oriented along the z-axis (this was chosen when it was created in Cubit)
waveguide = Waveguide3D.construct_simple("meshes/rectangular_waveguide_12000tets_example.inp", 4)

# Solve the FEM problem. This constructs the matrix equation and solves it by taking a matrix inverse
waveguide.solve()

# Determine how long the geometry is along the z-axis. z_max and z_min are attributes of the Waveguide3D object
y_length = waveguide.y_max - waveguide.y_min

# Produce some field results in the XY plane offset halfway along the geometry. This produces a matplotlib figure.
# On Bell, the interactive mode of matplotlib does not seem to work (it won't pop up an figure to interact with)
waveguide.plot_fields(plane="xz", offset=y_length/2, phase=0)

# matplotlib still works without interactive mode, but you will need to save the the results to view them
plt.savefig("fields_plot_simple.png")
