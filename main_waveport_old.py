import matplotlib.pyplot as plt
import numpy as np
from waveport.waveport import Waveguide3D
from scipy.constants import *
from math import atan2, floor
from waveport.util import plot_phases_csv, plot_csv

vn, vp, pn = ["TetrahedronsVacuum", "TetrahedronsSubstrate"], [1, 4.5], "PECWalls"
abcn = "ABC"
ipn, ipbn, ipp = ["InputPortVacuum", "InputPortSubstrate"], "InPortPEC", [1, 4.5]
opn, opbn, opp = ["OutputPortVacuum", "OutputPortSubstrate"], "OutPortPEC", [1, 4.5]
p1ip, p2ip = np.array([0, 0.0008]), np.array([0, 0])
p1op, p2op = np.array([0, 0.0008]), np.array([0, 0])
# Generate S21 results for a range of k0 (frequencies)
num_freqs = 7
microstrip = True
if microstrip:
    freqs = np.linspace(1E6, 100E6, num_freqs)
    k0s = freqs * 2 * pi / c
else:
    k0s = np.linspace(4, 4.5, num_freqs)
    freqs = k0s / 2 / pi * c
s21s = np.zeros([num_freqs])
s11s = np.zeros([num_freqs])
phases = np.zeros([num_freqs])
expected_phases = np.zeros([num_freqs])
for i, k0 in enumerate(k0s):
    # Fix a small issue that occurs at that particular frequency
    if i == 0 and microstrip:
        k0 += 0.0000000001
    print(f"Starting iteration {i+1} of {len(k0s)} for k0 = {k0}")
    if microstrip:
        waveguide = Waveguide3D("../cubit_meshes/fixed_pec_mesh.inp", k0, vn, vp, pn, abcn, ipn, ipbn, ipp, opn, opbn,
                                opp, p1ip, p2ip, p1op, p2op)
        index = 0
        while waveguide.input_port.get_selected_beta() > 10:
            index += 1
            waveguide.input_port.set_mode_index(index)
        index = 0
        while waveguide.output_port.get_selected_beta() > 10:
            index += 1
            waveguide.output_port.set_mode_index(index)
        waveguide.solve(index)
    else:
        waveguide = Waveguide3D.construct_simple(
            "meshes/rectangular_waveguide_12000tets_correct_orientation_20220630.inp", k0)
        waveguide.solve()
    lam = 2 * pi / waveguide.input_port.get_selected_beta()
    z_length = waveguide.z_max - waveguide.z_min
    expected_phases[i] = (z_length / lam * 360) % 360
    s11 = waveguide.compute_s11()
    s11s[i] = abs(s11)
    s21 = waveguide.compute_s21()
    s21s[i] = abs(s21)
    phases[i] = atan2(s21.imag, s21.real) / 2 / pi * 360
    if phases[i] < 0:
        phases[i] += 360
    print(f"Finished iteration {i+1} of {len(k0s)} for k0 = {k0}")
    print(f"S21 Phases: {phases}")
    print(f"S21 Mags: {s21s}")
    print(f"S11 Mags: {s11s}")

if microstrip:
    plt.figure()
    plt.plot(freqs/1E6, phases)
    plot_phases_csv("../scratch/microstrip_s21_phases.csv", 1, 1, True, False)
    plt.xlabel("Frequency (MHz)")
    plt.ylabel("Phase (degrees)")
    plt.xlim(-5, 105)
    plt.xticks([1] + [(i+1)*20 for i in range(5)])
    plt.ylim(350, 362)
    plt.legend(["FEM Code", "HFSS"])
    plt.tight_layout()
    plt.savefig("microstrip_s21_phases.png")
    plt.close()

    plt.figure()
    plt.plot(freqs/1E6, s21s)
    plot_csv("../scratch/microstrip_s21_mag.csv", 1, 1, True, False)
    plt.xlabel("Frequency (MHz)")
    plt.ylabel("Mag(S21)")
    plt.xlim(-5, 105)
    plt.xticks([1] + [(i+1)*20 for i in range(5)])
    plt.ylim(0.80, 1.02)
    plt.tight_layout()
    plt.legend(["FEM Code", "HFSS"])
    plt.savefig("microstrip_s21_mag.png")
    plt.close()

    plt.figure()
    plt.plot(freqs/1E6, s11s)
    plot_csv("../scratch/microstrip_s11_mag.csv", 1, 1, True, False)
    plt.xlabel("Frequency (MHz)")
    plt.ylabel("Mag(S11)")
    plt.xlim(-5, 105)
    plt.xticks([1] + [(i+1)*20 for i in range(5)])
    plt.ylim(0, 0.2)
    plt.tight_layout()
    plt.legend(["FEM Code", "HFSS"])
    plt.savefig("microstrip_s11_mag.png")
    plt.close()
else:
    # Make units in terms of MHz
    xlabel = "Frequency (MHz)"
    freqs /= 1E6
    plt.figure()
    plt.plot(freqs, s21s)
    plt.plot(freqs, np.loadtxt("../scratch/low_freq_s21s.csv"), "--")
    plt.xlabel(xlabel)
    plt.ylabel("Mag(S21)")
    plt.legend(["Field-Based", "Potential-Based"])
    # plt.xlim(-5, 105)
    # plt.xticks([1] + [(i + 1) * 20 for i in range(10)])
    plt.ylim(0.90, 1)
    plt.tight_layout()
    plt.savefig("field_vs_potential_s21.png")
    plt.close()

    plt.figure()
    plt.plot(freqs, s11s)
    plt.plot(freqs, np.loadtxt("../scratch/low_freq_s11s.csv"), "--")
    plt.xlabel(xlabel)
    plt.ylabel("Mag(S11)")
    plt.legend(["Field-Based", "Potential-Based"])
    # plt.xlim(-5, 105)
    # plt.xticks([1] + [(i + 1) * 20 for i in range(10)])
    plt.ylim(0, 0.1)
    plt.tight_layout()
    plt.savefig("field_vs_potential_s11.png")
    plt.close()

# # plt.plot(k0s, expected_phases)
# plt.savefig("s21_phases_final.png")
# plt.close()
# The most common test so far:
# waveguide = Waveguide3D("rectangular_waveguide_finer_20220625.inp", 4)
waveguide = Waveguide3D.construct_simple("meshes/rectangular_waveguide_12000tets_correct_orientation_20220630.inp", 4)
waveguide.solve()
# waveguide = Waveguide3D("rectangular_waveguide_12000tets_correct_orientation_rotated_20220705.inp", 8)
# Create the microstrip-line simulated at f = 1 MHz
# waveguide = Waveguide3D("microstrip_line_44000tets_20220710.inp", 0.0209584502195, vn, vp, pn, ipn, ipbn, ipp, opn, opbn, opp)
# waveguide = Waveguide3D("microstrip_line_44000tets_20220710.inp", 4, vn, vp, pn, None, ipn, ipbn, ipp, opn, opbn, opp, p1ip, p2ip, p1op, p2op)
# waveguide = Waveguide3D("microstrip_line_44000tets_20220710.inp", 220, vn, vp, pn, None, ipn, ipbn, ipp, opn, opbn, opp, p1ip, p2ip, p1op, p2op)
# waveguide = Waveguide3D("microstrip_line_44000tets_with_abc_20220711.inp", 220, vn, vp, pn, abcn, ipn, ipbn, ipp, opn, opbn, opp, p1ip, p2ip, p1op, p2op)
# waveguide = Waveguide3D("../cubit_meshes/fixed_abc_mesh.inp", 220, vn, vp, pn, abcn, ipn, ipbn, ipp, opn, opbn, opp, p1ip, p2ip, p1op, p2op)
# waveguide = Waveguide3D("../cubit_meshes/fixed_pec_mesh.inp", 2.09585, vn, vp, pn, abcn, ipn, ipbn, ipp, opn, opbn, opp, p1ip, p2ip, p1op, p2op)

# TYPICAL SIMULATION
# waveguide = Waveguide3D("../cubit_meshes/fixed_pec_mesh.inp", 0.02095845, vn, vp, pn, abcn, ipn, ipbn, ipp, opn, opbn, opp, p1ip, p2ip, p1op, p2op)

# waveguide = Waveguide3D("../cubit_meshes/fixed_pec_mesh.inp", 220, vn, vp, pn, abcn, ipn, ipbn, ipp, opn, opbn, opp, p1ip, p2ip, p1op, p2op)

s11 = waveguide.compute_s11()
s21 = waveguide.compute_s21()
print(f"S11: {abs(s11)} @ {atan2(s11.imag, s11.real) / 2 / pi * 360} deg")
print(f"S21: {abs(s21)} @ {atan2(s21.imag, s21.real) / 2 / pi * 360} deg")
z_length = waveguide.z_max - waveguide.z_min
y_length = waveguide.y_max - waveguide.y_min

waveguide.Ex = None
vmin = -0.03
vmax = 0.03
num_phases = 50
for i in range(num_phases):
    waveguide.plot_fields(plane="xy", offset=z_length/2, phase=i*2*pi/num_phases, use_cached_fields=True, vmin=vmin, vmax=vmax)
    # plt.savefig(f"images/te10_planexy_{floor(i/10)}{i%10}")
    plt.savefig(f"images/tem_planexy_{floor(i/10)}{i%10}")
    plt.close()
vmin = -0.12
vmax = 0.12
waveguide.Ex = None
for i in range(num_phases):
    waveguide.plot_fields(plane="xz", offset=y_length/2 + 0.0004, phase=i*2*pi/num_phases, vmin=vmin, vmax=vmax, use_cached_fields=True)
    plt.savefig(f"images/te10_planexz_y0p25_{floor(i/10)}{i%10}")
    plt.close()

waveguide.plot_one_field("Ex", offset=z_length/2, normalize=True)
plt.savefig("Ex_planexy_center.png")
plt.close()
waveguide.plot_one_field("Ey", offset=z_length/2, normalize=True)
plt.savefig("Ey_planexy_center.png")
plt.close()
print("Done creating plots")
