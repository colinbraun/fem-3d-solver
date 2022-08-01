import matplotlib.pyplot as plt
import numpy as np
from waveport.waveport import Waveguide3D
from waveport.waveport_lfbd import Waveguide3DLFBD
from scipy.constants import *
from math import atan2, floor

vn, vp, pn = ["TetrahedronsVacuum", "TetrahedronsSubstrate"], [1, 4.5], "PECWalls"
abcn = "ABC"
ipn, ipbn, ipp = ["InputPortVacuum", "InputPortSubstrate"], "InPortPEC", [1, 4.5]
opn, opbn, opp = ["OutputPortVacuum", "OutputPortSubstrate"], "OutPortPEC", [1, 4.5]
p1ip, p2ip = np.array([0, 0.0008]), np.array([0, 0])
p1op, p2op = np.array([0, 0.0008]), np.array([0, 0])
waveguide = Waveguide3DLFBD.construct_simple("rectangular_waveguide_12000tets_correct_orientation_20220630.inp", 4)
waveguide.solve()
z_length = waveguide.z_max - waveguide.z_min
y_length = waveguide.y_max - waveguide.y_min
waveguide.plot_fields(plane="xy", offset=z_length / 2, phase=0, use_cached_fields=True)
plt.savefig("lfbd_fields.png")
