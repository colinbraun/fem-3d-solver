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
# waveguide = Waveguide3DLFBD.construct_simple("rectangular_waveguide_12000tets_correct_orientation_20220630.inp", 4)
waveguide = Waveguide3D("meshes/microstrip_line_44000tets_pec_walls_example.inp", 2.095845, vn, vp, pn, abcn, ipn, ipbn, ipp, opn, opbn, opp, p1ip, p2ip, p1op, p2op)
waveguide.solve()
z_length = waveguide.z_max - waveguide.z_min
y_length = waveguide.y_max - waveguide.y_min
# waveguide.plot_fields(plane="xy", offset=z_length / 2, phase=0, use_cached_fields=True)
# plt.savefig("lfbd_fields.png")
waveguide.plot_one_field("Ex", offset=z_length/2, normalize=True, vmin=-0.6, vmax=1)
plt.savefig("Ex_planexy_center_matched_range.png")
plt.close()
waveguide.plot_one_field("Ey", offset=z_length/2, normalize=True, vmin=-0.4, vmax=1)
plt.savefig("Ey_planexy_center_matched_range.png")
plt.close()
print("Done creating plots")

# waveguide.Ex = None
# vmin = None
# vmax = None
# num_phases = 50
# for i in range(num_phases):
#     waveguide.plot_fields(plane="xy", offset=z_length/2, phase=i*2*pi/num_phases, use_cached_fields=True, vmin=vmin, vmax=vmax)
#     # plt.savefig(f"images/te10_planexy_{floor(i/10)}{i%10}")
#     plt.savefig(f"images/tem_planexy_{floor(i/10)}{i%10}")
#     plt.close()
# vmin = None
# vmax = None
# waveguide.Ex = None
# for i in range(num_phases):
#     waveguide.plot_fields(plane="xz", offset=y_length/2 + 0.0004, phase=i*2*pi/num_phases, vmin=vmin, vmax=vmax, use_cached_fields=True)
#     plt.savefig(f"images/te10_planexz_y0p25_{floor(i/10)}{i%10}")
#     plt.close()
