from waveport.util import load_mesh, Edge, quad_eval, quad_sample_points, where, plot_csv, plot_phases_csv
import numpy as np
from iwaveguide.waveguide import Waveguide
from scipy.linalg import inv
import scipy.sparse as sparse
from scipy.sparse.linalg import gmres, spilu
import matplotlib.pyplot as plt
from math import floor, e, pi, atan2, sqrt
from scipy.constants import c, mu_0, epsilon_0
import time

# Turn on interactive plotting
plt.ion()


class Waveguide3D:
    """Class representing a 3D waveguide with an input port and output port (not necessarily same size or shape)"""
    def __init__(self, filename, k0, volume_names, volume_permittivities, pec_name, abc_name, ip_surface_names, ip_boundary_name, ip_permittivities, op_surface_names, op_boundary_name, op_permittivities, p1_ip=None, p2_ip=None, p1_op=None, p2_op=None):
        """
        Constructor of Waveguide3D.
        :param filename: A string giving the path to the .inp file to be loaded.
        :param k0: The k0 value (frequency) to operate at.
        :param volume_names: An iterable containing the names of the blocks from cubit for each volume. This allows for
        different parts of the geometry to have different permittivity values.
        :param volume_permittivities: A list of the relative permittivities of the material for the volumes.
        :param pec_name: A string of the name of the block from cubit containing the PEC surface information.
        If ``None`` is passed, it will be assumed that there are no PEC walls.
        :param abc_name: A string of the name of the block from cubit containing the Absorbing Boundary surface information.
        If ``None`` is passed, it will be assumed that there are no absorbing boundary conditions.
        :param ip_surface_names: A list containing the names of each surface of the mesh in the input port.
        :param ip_boundary_name: A string of the name of the boundary nodes of the mesh for the input port.
        :param ip_permittivities: A list of the relative permittivities of the material for the input port surfaces.
        :param op_surface_names: A list containing the names of each surface of the mesh in the output port.
        :param op_boundary_name: A string of the name of the boundary nodes of the mesh for the output port.
        :param op_permittivities: A list of the relative permittivities of the material for the output port surfaces.
        :param p1_ip: If an integration line is desired, ``p1_ip`` is the starting point of the line integral on the
        input port. Default = None.
        :param p2_ip: If an integration line is desired, ``p2_ip`` is the starting point of the line integral on the
        input port. Default = None.
        :param p1_op: If an integration line is desired, ``p1_op`` is the starting point of the line integral on the
        output port. Default = None.
        :param p2_op: If an integration line is desired, ``p2_op`` is the starting point of the line integral on the
        output port. Default = None.
        """
        # Load the mesh
        start_time = time.time()
        self.all_nodes, self.tetrahedrons, self.tets_node_ids, self.all_edges, self.boundary_pec_edge_numbers, self.abc_edge_numbers, self.boundary_input_edge_numbers, self.boundary_output_edge_numbers, self.remap_edge_nums, self.all_edges_map, self.boundary_input_tets, self.boundary_output_tets = load_mesh(filename, volume_names, volume_permittivities, pec_name, abc_name, ip_surface_names, ip_boundary_name, ip_permittivities, op_surface_names, op_boundary_name, op_permittivities)
        print(f"Mesh loaded in {time.time() - start_time} seconds")
        # Find the minimum and maximum x, y, and z values of the waveguide.
        self.x_min = np.amin(self.all_nodes[:, 0])
        self.x_max = np.amax(self.all_nodes[:, 0])
        self.y_min = np.amin(self.all_nodes[:, 1])
        self.y_max = np.amax(self.all_nodes[:, 1])
        self.z_min = np.amin(self.all_nodes[:, 2])
        self.z_max = np.amax(self.all_nodes[:, 2])
        self.k0 = k0
        self.f = k0 / 2 / pi * c
        # The angular frequency of operation (from k0)
        self.omega = k0 / sqrt(mu_0 * epsilon_0)
        # Construct empty K and b matrices
        self.K = np.zeros([len(self.remap_edge_nums), len(self.remap_edge_nums)], dtype=complex)
        self.b = np.zeros([len(self.remap_edge_nums)], dtype=complex)
        # Create a Waveguide object of the input port
        self.input_port = Waveguide(filename, ip_surface_names, ip_boundary_name, ip_permittivities, p1_ip, p2_ip)
        # Set its mode to be the specified propagating mode (0 -> TE10, 4 -> TM11)
        self.input_port.set_mode_index(0)
        # Solve the waveguide for k0 = 4 (this will work for TE10 mode excitation, but won't excite other modes)
        self.input_port.solve_k0(k0)
        # TODO: Change this to have an output port that differs from the input port. Will have to solve it too.
        # self.output_port = self.input_port
        self.output_port = Waveguide(filename, op_surface_names, op_boundary_name, op_permittivities, p1_op, p2_op)
        self.output_port.set_mode_index(0)
        self.output_port.solve_k0(k0)
        # Create an empty coefficient array. This is loaded by solve()
        self.edge_coefficients = np.zeros([len(self.remap_edge_nums), len(self.remap_edge_nums)])
        # Placeholders
        self.Ex, self.Ey, self.Ez = None, None, None

    @staticmethod
    def construct_simple(filename, k0=4., permittivity=1., p1_ip=None, p2_ip=None, p1_op=None, p2_op=None):
        """
        Simple construction of a ``Waveguide3D`` object, avoiding needing to specify the exhausting details.
        This is essentially the old constructor.
        :param filename: A string giving the path to the .inp file to be loaded.
        :param k0: The k0 value (frequency) to operate at.
        :param permittivity: The relative permittivity of the material
        :param p1_ip: If an integration line is desired, ``p1_ip`` is the starting point of the line integral on the
        input port. Default = None.
        :param p2_ip: If an integration line is desired, ``p2_ip`` is the starting point of the line integral on the
        input port. Default = None.
        :param p1_op: If an integration line is desired, ``p1_op`` is the starting point of the line integral on the
        output port. Default = None.
        :param p2_op: If an integration line is desired, ``p2_op`` is the starting point of the line integral on the
        output port. Default = None.
        """
        return Waveguide3D(filename, k0, ["Tetrahedrons"], [permittivity], "PECWalls", None, ["InputPort"], "InPortBoundary", [permittivity], ["OutputPort"], "OutPortBoundary", [permittivity], p1_ip, p2_ip, p1_op, p2_op)

    def solve(self, mode=0):
        """
        Solve for the fields in the Waveguide3D object.
        :param mode: The mode to excite at the input port. Default = 0 (lowest cut-off frequency mode).
        :return: The edge coefficients to be applied to the interpolating functions (indexed by global edge number).
        """
        start_time = time.time()
        print("Begin constructing equation matrix")
        # Create the b vector
        self.b = self.generate_b_vector(self.input_port, mode)
        # Iterate over the tetrahedrons and construct the K and b matrices (this concept comes from Jin page 454)
        for tet in self.tetrahedrons:

            # Compute the x-mean, y-mean, and z-mean of the points that make up the tetrahedron (needed later)
            x_mean = np.average(tet.points[:, 0])
            y_mean = np.average(tet.points[:, 1])
            z_mean = np.average(tet.points[:, 2])

            # Iterate over the edges of the tetrahedron
            for edgei in tet.edges:
                # Skip over PEC walls
                if edgei in self.boundary_pec_edge_numbers:
                    continue
                # Get a hold of the Edge object
                edge1 = self.all_edges[edgei]

                # The simplex constants for nodes i and j of edge l
                indices_l = [np.argwhere(tet.nodes == edge1.node1)[0][0], np.argwhere(tet.nodes == edge1.node2)[0][0]]
                # Necessary constants from NASA paper eqs. 163-172
                a_il, a_jl = tet.simplex_consts[indices_l]
                Axl = a_il[0]*a_jl[1] - a_il[1]*a_jl[0]
                Bxl = a_il[2]*a_jl[1] - a_il[1]*a_jl[2]
                Cxl = a_il[3]*a_jl[1] - a_il[1]*a_jl[3]
                Ayl = a_il[0]*a_jl[2] - a_il[2]*a_jl[0]
                Byl = a_il[1]*a_jl[2] - a_il[2]*a_jl[1]
                Cyl = a_il[3]*a_jl[2] - a_il[2]*a_jl[3]
                Azl = a_il[0]*a_jl[3] - a_il[3]*a_jl[0]
                Bzl = a_il[1]*a_jl[3] - a_il[3]*a_jl[1]
                Czl = a_il[2]*a_jl[3] - a_il[3]*a_jl[2]

                # Iterate over the edges of the tetrahedron
                for edgej in tet.edges:
                    # Skip over PEC walls
                    if edgej in self.boundary_pec_edge_numbers:
                        continue
                    edge2 = self.all_edges[edgej]
                    node_ik, node_jk = self.all_nodes[edge2.node1], self.all_nodes[edge2.node2]

                    # Find the indices of the edge of interest
                    indices_k = [np.argwhere(tet.nodes == edge2.node1)[0][0], np.argwhere(tet.nodes == edge2.node2)[0][0]]
                    # The simplex constants for nodes i and j of edge l
                    a_ik, a_jk = tet.simplex_consts[indices_k]
                    # Necessary constants from NASA paper eqs. 163-172
                    Axk = a_ik[0] * a_jk[1] - a_ik[1] * a_jk[0]
                    Bxk = a_ik[2] * a_jk[1] - a_ik[1] * a_jk[2]
                    Cxk = a_ik[3] * a_jk[1] - a_ik[1] * a_jk[3]
                    Ayk = a_ik[0] * a_jk[2] - a_ik[2] * a_jk[0]
                    Byk = a_ik[1] * a_jk[2] - a_ik[2] * a_jk[1]
                    Cyk = a_ik[3] * a_jk[2] - a_ik[2] * a_jk[3]
                    Azk = a_ik[0] * a_jk[3] - a_ik[3] * a_jk[0]
                    Bzk = a_ik[1] * a_jk[3] - a_ik[3] * a_jk[1]
                    Czk = a_ik[2] * a_jk[3] - a_ik[3] * a_jk[2]
                    # If working with input port edges, need to do different integral
                    if edgei in self.boundary_input_edge_numbers and edgej in self.boundary_input_edge_numbers:
                        # Perform the (n_hat x N_i) \dot (n_hat x N_j) integral
                        # For n_hat = z_hat, this turns out to be the same as N_i \dot N_j, so integral is same as inhomo wg

                        # Need the third point that makes up the face of the input port (first 2 are in the edge)
                        found_edge_no = -1
                        # Search for a different edge on the input port from this tetrahedral element
                        for edge in tet.edges:
                            if edge == edgej:
                                continue
                            elif edge in self.boundary_input_edge_numbers:
                                # We found an edge containing the third node, note it
                                found_edge_no = edge
                                break
                        # It is possible we do not find another edge that lies on the input port. We are checking for this here.
                        if found_edge_no != -1:
                            found_edge = self.all_edges[found_edge_no]
                            nodes = np.unique(
                                np.array([node_ik, node_jk, self.all_nodes[found_edge.node1], self.all_nodes[found_edge.node2]]), axis=0)
                            if len(nodes) != 3:
                                raise RuntimeError("Did not find 3 nodes on surface triangle")
                            # Get two vectors that lie in the triangular surface
                            v1 = nodes[1] - nodes[0]
                            v2 = nodes[2] - nodes[0]
                            # We SHOULD find the n_hat that points out of the surface, but inward pointing will cancel anyway
                            n_hat = np.cross(v1, v2)
                            # Normalize n_hat
                            n_hat /= np.linalg.norm(n_hat)
                            sample_points = quad_sample_points(3, nodes[0], nodes[1], nodes[2])
                            N_i = np.zeros([len(sample_points), 3])
                            N_i[:, 0] = Axl + Bxl * sample_points[:, 1] + Cxl * sample_points[:, 2]
                            N_i[:, 1] = Ayl + Byl * sample_points[:, 0] + Cyl * sample_points[:, 2]
                            N_i[:, 2] = Azl + Bzl * sample_points[:, 0] + Czl * sample_points[:, 1]
                            N_i = N_i * edge1.length / 36 / tet.volume**2
                            n_hat_x_ni = np.array([np.cross(n_hat, vec) for vec in N_i])
                            N_j = np.zeros([len(sample_points), 3])
                            N_j[:, 0] = Axk + Bxk * sample_points[:, 1] + Cxk * sample_points[:, 2]
                            N_j[:, 1] = Ayk + Byk * sample_points[:, 0] + Cyk * sample_points[:, 2]
                            N_j[:, 2] = Azk + Bzk * sample_points[:, 0] + Czk * sample_points[:, 1]
                            N_j = N_j * edge2.length / 36 / tet.volume**2
                            n_hat_x_nj = np.array([np.cross(n_hat, vec) for vec in N_j])
                            dotted = n_hat_x_ni[:, 0] * n_hat_x_nj[:, 0] + n_hat_x_ni[:, 1] * n_hat_x_nj[:, 1] + n_hat_x_ni[:, 2] * n_hat_x_nj[:, 2]
                            ni_dot_nj = np.reshape(N_i[:, 0] * N_j[:, 0] + N_i[:, 1] * N_j[:, 1] + N_i[:, 2] * N_j[:, 2], [len(sample_points), 1])
                            dotted = np.reshape(dotted, [len(dotted), 1])
                            integral = quad_eval(nodes[0], nodes[1], nodes[2], dotted)
                            self.K[self.remap_edge_nums[edgei], self.remap_edge_nums[edgej]] += integral * self.input_port.get_selected_beta() * 1j

                    # If working with output port edges, need to do different integral
                    elif edgei in self.boundary_output_edge_numbers and edgej in self.boundary_output_edge_numbers:
                        # Perform the (n_hat x N_i) \dot (n_hat x N_j) integral
                        # For n_hat = z_hat, this turns out to be the same as N_i \dot N_j, so integral is same as inhomo wg

                        # Need the third point that makes up the face of the input port (first 2 are in the edge)
                        found_edge_no = -1
                        # Search for a different edge on the input port from this tetrahedral element
                        for edge in tet.edges:
                            if edge == edgej:
                                continue
                            elif edge in self.boundary_output_edge_numbers:
                                # We found an edge containing the third node, note it
                                found_edge_no = edge
                                break
                        # It is possible we do not find another edge that lies on the input port. We are checking for this here.
                        if found_edge_no != -1:
                            found_edge = self.all_edges[found_edge_no]
                            nodes = np.unique(
                                np.array([node_ik, node_jk, self.all_nodes[found_edge.node1], self.all_nodes[found_edge.node2]]), axis=0)
                            if len(nodes) != 3:
                                raise RuntimeError("Did not find 3 nodes on surface triangle")
                            # Get two vectors that lie in the triangular surface
                            v1 = nodes[1] - nodes[0]
                            v2 = nodes[2] - nodes[0]
                            # We SHOULD find the n_hat that points out of the surface, but inward pointing will cancel anyway
                            n_hat = np.cross(v1, v2)
                            # Normalize n_hat
                            n_hat /= np.linalg.norm(n_hat)
                            sample_points = quad_sample_points(3, nodes[0], nodes[1], nodes[2])
                            N_i = np.zeros([len(sample_points), 3])
                            N_i[:, 0] = Axl + Bxl * sample_points[:, 1] + Cxl * sample_points[:, 2]
                            N_i[:, 1] = Ayl + Byl * sample_points[:, 0] + Cyl * sample_points[:, 2]
                            N_i[:, 2] = Azl + Bzl * sample_points[:, 0] + Czl * sample_points[:, 1]
                            N_i = N_i * edge1.length / 36 / tet.volume**2
                            n_hat_x_ni = np.array([np.cross(n_hat, vec) for vec in N_i])
                            N_j = np.zeros([len(sample_points), 3])
                            N_j[:, 0] = Axk + Bxk * sample_points[:, 1] + Cxk * sample_points[:, 2]
                            N_j[:, 1] = Ayk + Byk * sample_points[:, 0] + Cyk * sample_points[:, 2]
                            N_j[:, 2] = Azk + Bzk * sample_points[:, 0] + Czk * sample_points[:, 1]
                            N_j = N_j * edge2.length / 36 / tet.volume**2
                            n_hat_x_nj = np.array([np.cross(n_hat, vec) for vec in N_j])
                            dotted = n_hat_x_ni[:, 0] * n_hat_x_nj[:, 0] + n_hat_x_ni[:, 1] * n_hat_x_nj[:, 1] + n_hat_x_ni[:, 2] * n_hat_x_nj[:, 2]
                            dotted = np.reshape(dotted, [len(dotted), 1])
                            integral = quad_eval(nodes[0], nodes[1], nodes[2], dotted)
                            self.K[self.remap_edge_nums[edgei], self.remap_edge_nums[edgej]] += self.output_port.get_selected_beta() * integral * 1j

                    elif edgei in self.abc_edge_numbers and edgej in self.abc_edge_numbers:
                        # Perform the (n_hat x N_i) \dot (n_hat x N_j) integral
                        # For n_hat = z_hat, this turns out to be the same as N_i \dot N_j, so integral is same as inhomo wg

                        # Need the third point that makes up the face of the input port (first 2 are in the edge)
                        found_edge_no = -1
                        # Search for a different edge on the input port from this tetrahedral element
                        for edge in tet.edges:
                            if edge == edgej:
                                continue
                            elif edge in self.boundary_output_edge_numbers:
                                # We found an edge containing the third node, note it
                                found_edge_no = edge
                                break
                        # It is possible we do not find another edge that lies on the input port. We are checking for this here.
                        if found_edge_no != -1:
                            found_edge = self.all_edges[found_edge_no]
                            nodes = np.unique(
                                np.array([node_ik, node_jk, self.all_nodes[found_edge.node1], self.all_nodes[found_edge.node2]]), axis=0)
                            if len(nodes) != 3:
                                raise RuntimeError("Did not find 3 nodes on surface triangle")
                            # Get two vectors that lie in the triangular surface
                            v1 = nodes[1] - nodes[0]
                            v2 = nodes[2] - nodes[0]
                            # We SHOULD find the n_hat that points out of the surface, but inward pointing will cancel anyway
                            n_hat = np.cross(v1, v2)
                            # Normalize n_hat
                            n_hat /= np.linalg.norm(n_hat)
                            sample_points = quad_sample_points(3, nodes[0], nodes[1], nodes[2])
                            N_i = np.zeros([len(sample_points), 3])
                            N_i[:, 0] = Axl + Bxl * sample_points[:, 1] + Cxl * sample_points[:, 2]
                            N_i[:, 1] = Ayl + Byl * sample_points[:, 0] + Cyl * sample_points[:, 2]
                            N_i[:, 2] = Azl + Bzl * sample_points[:, 0] + Czl * sample_points[:, 1]
                            N_i = N_i * edge1.length / 36 / tet.volume**2
                            n_hat_x_ni = np.array([np.cross(n_hat, vec) for vec in N_i])
                            N_j = np.zeros([len(sample_points), 3])
                            N_j[:, 0] = Axk + Bxk * sample_points[:, 1] + Cxk * sample_points[:, 2]
                            N_j[:, 1] = Ayk + Byk * sample_points[:, 0] + Cyk * sample_points[:, 2]
                            N_j[:, 2] = Azk + Bzk * sample_points[:, 0] + Czk * sample_points[:, 1]
                            N_j = N_j * edge2.length / 36 / tet.volume**2
                            n_hat_x_nj = np.array([np.cross(n_hat, vec) for vec in N_j])
                            dotted = n_hat_x_ni[:, 0] * n_hat_x_nj[:, 0] + n_hat_x_ni[:, 1] * n_hat_x_nj[:, 1] + n_hat_x_ni[:, 2] * n_hat_x_nj[:, 2]
                            dotted = np.reshape(dotted, [len(dotted), 1])
                            integral = quad_eval(nodes[0], nodes[1], nodes[2], dotted)
                            self.K[self.remap_edge_nums[edgei], self.remap_edge_nums[edgej]] += 2 * pi * self.f / sqrt(epsilon_0 * mu_0) * integral * 1j
                            continue

                    # Compute the curl(N_i) \dot curl(N_j) part of K_ij
                    curl_dot_curl_part = edge1.length*edge2.length / 324 / tet.volume**3 * (Czl*Czk + Cxl*Cxk + Byl*Byk)

                    # Build constants for N_i \dot N_j integral (part of K_ij and b_i)
                    # These constants come from NASA paper eq. 182
                    I1 = Axl*Axk + Ayl*Ayk + Azl*Azk
                    I2 = (Ayl*Byk + Ayk*Byl + Azl*Bzk + Azk*Bzl) * x_mean
                    I3 = (Axl*Bxk + Axk*Bxl + Azl*Czk + Azk*Czl) * y_mean
                    I4 = (Axl*Cxk + Axk*Cxl + Ayl*Cyk + Ayk*Cyl) * z_mean
                    I5 = 1/20 * (Bzl*Czk + Bzk*Czl) * (sum(xi*yi for xi, yi, zi in tet.points) + 16*x_mean*y_mean)
                    I6 = 1/20 * (Bxl*Cxk + Bxk*Cxl) * (sum(yi*zi for xi, yi, zi in tet.points) + 16*y_mean*z_mean)
                    I7 = 1/20 * (Byl*Cyk + Byk*Cyl) * (sum(xi*zi for xi, yi, zi in tet.points) + 16*x_mean*z_mean)
                    I8 = 1/20 * (Byl*Byk + Bzl*Bzk) * (sum(xi*xi for xi, yi, zi in tet.points) + 16*x_mean*x_mean)
                    I9 = 1/20 * (Bxl*Bxk + Czl*Czk) * (sum(yi*yi for xi, yi, zi in tet.points) + 16*y_mean*y_mean)
                    I10 = 1/20 * (Cxl*Cxk + Cyl*Cyk) * (sum(zi*zi for xi, yi, zi in tet.points) + 16*z_mean*z_mean)
                    i_sum = I1 + I2 + I3 + I4 + I5 + I6 + I7 + I8 + I9 + I10
                    dot_part = -self.input_port.k0**2 * tet.permittivity * edge1.length * edge2.length / 1296 / tet.volume**3 * i_sum
                    k_value = curl_dot_curl_part + dot_part
                    index1 = self.remap_edge_nums[edgei]
                    index2 = self.remap_edge_nums[edgej]
                    self.K[index1, index2] += k_value
        print(f"Constructed equation matrix in {time.time() - start_time} seconds")
        print("Solving equation matrix")
        start_time = time.time()
        # self.s_A = sparse.csr_matrix(self.K)
        # precon = spilu(self.s_A)
        # self.edge_coefficients, info = gmres(self.s_A, self.b, tol=1E-6, callback=gmres_cb, callback_type="pr_norm", M=precon)
        # print(f"gmres info: {info}")
        self.edge_coefficients = np.dot(inv(self.K), self.b)
        print(f"Solved equation matrix in {time.time() - start_time} seconds")
        return self.edge_coefficients

    def compute_s11(self):
        """
        Integrate the input port, generating the S11 value. The integrals follow the Comsol webpage
        :return: The evaluated S-parameter.
        """
        mode = 0
        # The numerator of the S11 parameter (will need to normalize this manually)
        integral1 = 0
        # The denominator of the S11 parameter (comes normalized be default)
        integral2 = self.integrate_port_profile(self.input_port, mode)

        for tet in self.boundary_input_tets:
            # Want to collect 2 edges that lie on the surface to find the 3 nodes that make it up
            found_edge_nos = []
            for edge in tet.edges:
                if edge in self.boundary_input_edge_numbers:
                    # We found an edge containing the third node, note it
                    found_edge_nos.append(edge)
                # Once 2 edges have been found, stop searching
                if len(found_edge_nos) == 2:
                    break
            # This should always be False, but check just in case.
            if len(found_edge_nos) != 2:
                raise RuntimeError("Did not find boundary face of boundary tetrahedron")

            found_edge1 = self.all_edges[found_edge_nos[0]]
            found_edge2 = self.all_edges[found_edge_nos[1]]
            # Get the 3 nodes that make up the triangle on the surface
            nodes = np.unique(np.array([self.all_nodes[found_edge1.node1], self.all_nodes[found_edge1.node2],
                                        self.all_nodes[found_edge2.node1], self.all_nodes[found_edge2.node2]]),
                              axis=0)
            # Ensure 3 were found
            if len(nodes) != 3:
                raise RuntimeError("Did not find 3 nodes on surface triangle")
            # -------------Perform quadrature on the 3 nodes that make up the triangle on the surface-----------
            # Generate the points to sample at for the quadrature
            sample_points = quad_sample_points(3, nodes[0], nodes[1], nodes[2])
            # Get the edge coefficients for each of the edges that make up this tetrahedron
            phis = [self.edge_coefficients[self.remap_edge_nums[edge_no]] if edge_no in self.remap_edge_nums else 0 for
                    edge_no in tet.edges]
            # Interpolate the field measured at the input port
            E_ip = np.array([np.array(tet.interpolate(phis, sample_point)) for sample_point in sample_points])
            # Get the incident field at the input port
            E_inc = np.array([self.input_port.get_field_at(p[0], p[1]) for p in sample_points])
            # Ec - E1 for use in the integral
            E_ip_m_E_inc = E_ip - E_inc
            # Get the complex conjugate of the incident field
            E_inc_conj = np.conjugate(E_inc)
            # Compute the dot product at each point
            values = np.reshape(E_ip_m_E_inc[:, 0] * E_inc_conj[:, 0] + E_ip_m_E_inc[:, 1] * E_inc_conj[:, 1] + E_ip_m_E_inc[:, 2] * E_inc_conj[:, 2],
                                 [len(sample_points), 1])
            # Compute the contribution to the integral in the numerator using quadrature
            integral1 += quad_eval(nodes[0], nodes[1], nodes[2], values)
        return integral1 / integral2

    def compute_s21(self):
        """
        Compute the S21 value.
        :return: The S21 value.
        """
        # The numerator of the result
        integral1 = 0
        mode = 0
        # The denominator of the result
        integral2 = self.integrate_port_profile(self.input_port, mode)

        # Iterate over the output port boundary tetrahedrons to calculate integral1
        for tet in self.boundary_output_tets:
            # Want to collect 2 edges that lie on the surface to find the 3 nodes that make it up
            found_edge_nos = []
            for edge in tet.edges:
                if edge in self.boundary_output_edge_numbers:
                    # We found an edge containing the third node, note it
                    found_edge_nos.append(edge)
                # Once 2 edges have been found, stop searching
                if len(found_edge_nos) == 2:
                    break
            # This should always be False, but check just in case.
            if len(found_edge_nos) != 2:
                raise RuntimeError("Did not find boundary face of boundary tetrahedron")

            found_edge1 = self.all_edges[found_edge_nos[0]]
            found_edge2 = self.all_edges[found_edge_nos[1]]
            # Get the 3 nodes that make up the triangle on the surface
            nodes = np.unique(np.array([self.all_nodes[found_edge1.node1], self.all_nodes[found_edge1.node2],
                                        self.all_nodes[found_edge2.node1], self.all_nodes[found_edge2.node2]]),
                              axis=0)
            # Ensure 3 were found
            if len(nodes) != 3:
                raise RuntimeError("Did not find 3 nodes on surface triangle")
            # -------------Perform quadrature on the 3 nodes that make up the triangle on the surface-----------
            # Generate the points to sample at for the quadrature
            sample_points = quad_sample_points(3, nodes[0], nodes[1], nodes[2])
            # Get the edge coefficients for each of the edges that make up this tetrahedron
            phis = [self.edge_coefficients[self.remap_edge_nums[edge]] if edge in self.remap_edge_nums else 0 for
                    edge in tet.edges]
            # Interpolate the field measured at the output port
            E2 = np.array([np.array(tet.interpolate(phis, sample_point)) for sample_point in sample_points])
            # Get the incident field at the output port
            E_inc = np.array([self.output_port.get_field_at(p[0], p[1]) for p in sample_points])
            E_inc_conj = np.conjugate(E_inc)
            # Compute the dot product at each point
            values1 = np.reshape(E2[:, 0] * E_inc_conj[:, 0] + E2[:, 1] * E_inc_conj[:, 1] + E2[:, 2] * E_inc_conj[:, 2],
                                 [len(sample_points), 1])
            # Compute the contribution to the integral in the numerator using quadrature
            integral1 += quad_eval(nodes[0], nodes[1], nodes[2], values1)
        return integral1 / integral2

    def generate_b_vector(self, port, mode):
        """
        Generate the b vector with ``port`` as the excited port for mode ``mode``.
        :param port: The iwaveguide.Waveguide object to generate the b vector for.
        :param mode: The mode index to excite the port with (0 is the lowest cut-off frequency mode).
        :return: The generated b vector (as a numpy array).
        """
        if port is self.input_port:
            boundary_tets = self.boundary_input_tets
            boundary_edge_numbers = self.boundary_input_edge_numbers
            # TODO: Make this identified properly instead of assuming
            z = self.z_min
        elif port is self.output_port:
            boundary_tets = self.boundary_output_tets
            boundary_edge_numbers = self.boundary_output_edge_numbers
            z = self.z_max
        else:
            raise ValueError("'port' parameter did not match the input or output port.")
        # Set the port into the correct mode
        port.set_mode_index(mode)
        # Create the empty b vector
        b = np.zeros([len(self.remap_edge_nums)], dtype=complex)

        # Iterate over each of the tetrahedrons
        for tet in boundary_tets:
            # We only care about tetrahedrons that lie on the port of interest
            # if tet not in boundary_tets:
            #     continue
            # Iterate over all the edges of t
            for edge_no in tet.edges:
                # Skip over PEC walls
                if edge_no in self.boundary_pec_edge_numbers:
                    continue
                # TODO: Think about checking if this edge number is in the boundary edge numbers
                # NOTE: The results of uncommenting the below 1 lines of code do not seem to change any results.
                # if edge_no not in boundary_edge_numbers:
                #     continue
                # Get a hold of the Edge object
                edge = self.all_edges[edge_no]
                # Get the indices of the tet.nodes list that correspond to the global node numbers of the edges
                indices_l = [np.argwhere(tet.nodes == edge.node1)[0][0], np.argwhere(tet.nodes == edge.node2)[0][0]]
                # The simplex constants for nodes i and j of edge l
                a_il, a_jl = tet.simplex_consts[indices_l]
                # Necessary constants from NASA paper eqs. 163-172
                Axl = a_il[0]*a_jl[1] - a_il[1]*a_jl[0]
                Bxl = a_il[2]*a_jl[1] - a_il[1]*a_jl[2]
                Cxl = a_il[3]*a_jl[1] - a_il[1]*a_jl[3]
                Ayl = a_il[0]*a_jl[2] - a_il[2]*a_jl[0]
                Byl = a_il[1]*a_jl[2] - a_il[2]*a_jl[1]
                Cyl = a_il[3]*a_jl[2] - a_il[2]*a_jl[3]
                Azl = a_il[0]*a_jl[3] - a_il[3]*a_jl[0]
                Bzl = a_il[1]*a_jl[3] - a_il[3]*a_jl[1]
                Czl = a_il[2]*a_jl[3] - a_il[3]*a_jl[2]
                # Need to find the third point that makes up the face of the input port
                found_edge_nos = []
                # Search for a different edge on the input port from this tetrahedral element
                for edge_no2 in tet.edges:
                    if edge_no2 in boundary_edge_numbers:
                        # We found an edge containing the third node, note it
                        found_edge_nos.append(edge_no2)
                    if len(found_edge_nos) == 2:
                        break
                # It is possible we do not find another edge that lies on the input port. We are checking for this here.
                if len(found_edge_nos) != 2:
                    raise RuntimeError("Did not find 2 surface edges for tetrahedron marked as having a face on edge")
                found_edge1 = self.all_edges[found_edge_nos[0]]
                found_edge2 = self.all_edges[found_edge_nos[1]]
                nodes = np.unique(np.array([self.all_nodes[found_edge1.node1], self.all_nodes[found_edge1.node2],
                                            self.all_nodes[found_edge2.node1], self.all_nodes[found_edge2.node2]]),
                                  axis=0)
                if len(nodes) != 3:
                    raise RuntimeError("Did not find 3 nodes on surface triangle")
                # Perform the N_i \dot K_N integral for the surface using gaussian quadrature
                sample_points = quad_sample_points(3, nodes[0], nodes[1], nodes[2])
                N_i = np.zeros([len(sample_points), 3])
                N_i[:, 0] = Axl + Bxl * sample_points[:, 1] + Cxl * sample_points[:, 2]
                N_i[:, 1] = Ayl + Byl * sample_points[:, 0] + Cyl * sample_points[:, 2]
                N_i[:, 2] = Azl + Bzl * sample_points[:, 0] + Czl * sample_points[:, 1]
                N_i = N_i * edge.length / 36 / tet.volume ** 2
                # Compute the x component of the edge interpolating function for each of the sample points
                # TODO: Consider if I should multiply by e**(-1j*beta*z)
                # Get the E_inc field at each of the sample points
                # E_inc = -2j * port.get_selected_beta() * np.array(
                #     [np.array(port.get_field_at(sample_point[0], sample_point[1])) for sample_point in
                #      sample_points])
                E_inc = -2j * port.get_selected_beta() * np.array(
                    [np.array(port.get_field_at(sample_point[0], sample_point[1])) for sample_point in
                     sample_points])
                # Use analytical formulas
                # a, b_ = self.x_max - self.x_min, self.y_max - self.y_min
                # kz = analytical.rect_wg_beta(a, b_, self.k0*c)
                # E_inc = -2j*kz*e**(-1j*kz*z) * np.array([analytical.rect_wg_field_at(p[0]-self.x_min, p[1]-self.y_min, a, b_, self.k0*c) for p in sample_points])
                # Compute the dot product at each point
                values = np.reshape(N_i[:, 0] * E_inc[:, 0] + N_i[:, 1] * E_inc[:, 1] + N_i[:, 2] * E_inc[:, 2],
                                    [len(sample_points), 1])
                integral = quad_eval(nodes[0], nodes[1], nodes[2], values)
                b[self.remap_edge_nums[edge_no]] += -integral
        return b

    def integrate_port_profile(self, port, mode=0):
        """
        Integrate the field profile of port mode index ``mode`` on port ``port``. This is done by integrating the
        incident field profile using quadrature.
        :param port: A iwaveguide.Waveguide object corresponding to the port to integrate.
        :param mode: The mode of the profile to integrate. Default = 0 (lowest cut-off frequency mode). W.I.P.
        :return: The result as a complex number (in most cases the imaginary part will be 0).
        """
        if port is self.input_port:
            boundary_tets = self.boundary_input_tets
            boundary_edge_numbers = self.boundary_input_edge_numbers
        elif port is self.output_port:
            boundary_tets = self.boundary_output_tets
            boundary_edge_numbers = self.boundary_output_edge_numbers
        else:
            raise ValueError("'port' parameter did not match the input or output port.")
        integral = 0
        # Iterate over the input port boundary tetrahedrons to generate integral2
        for tet in boundary_tets:
            # Want to collect 2 edges that lie on the surface to find the 3 nodes that make it up
            found_edge_nos = []
            for edge in tet.edges:
                if edge in boundary_edge_numbers:
                    # We found an edge containing the third node, note it
                    found_edge_nos.append(edge)
                # Once 2 edges have been found, stop searching
                if len(found_edge_nos) == 2:
                    break
            # This should always be True, but check just in case.
            if len(found_edge_nos) != 2:
                raise RuntimeError("Did not find boundary face of boundary tetrahedron")
            found_edge1 = self.all_edges[found_edge_nos[0]]
            found_edge2 = self.all_edges[found_edge_nos[1]]
            # Get the 3 nodes that make up the triangle on the surface
            nodes = np.unique(np.array([self.all_nodes[found_edge1.node1], self.all_nodes[found_edge1.node2],
                                        self.all_nodes[found_edge2.node1], self.all_nodes[found_edge2.node2]]),
                              axis=0)
            # Ensure 3 were found
            if len(nodes) != 3:
                raise RuntimeError("Did not find 3 nodes on surface triangle")
            # -------------Perform quadrature on the 3 nodes that make up the triangle on the surface-----------
            # Generate the points to sample at for the quadrature
            sample_points = quad_sample_points(3, nodes[0], nodes[1], nodes[2])
            # Calculate the incident field at the port
            E_inc = np.array([port.get_field_at(p[0], p[1]) for p in sample_points])
            E_inc_conj = np.conjugate(E_inc)
            # Compute the dot product at each point
            values = np.reshape(E_inc[:, 0] * E_inc_conj[:, 0] + E_inc[:, 1] * E_inc_conj[:, 1] + E_inc[:, 2] * E_inc_conj[:, 2],
                                 [len(sample_points), 1])
            # Compute the integral in the denominator using quadrature
            integral += quad_eval(nodes[0], nodes[1], nodes[2], values)
        return integral

    def get_fields_in_plane(self, num_axis1_points=100, num_axis2_points=100, plane="xy", offset=0.1):
        """
        Compute the fields in the specified plane.
        :param num_axis1_points: The number of points to compute the fields for along the first axis in the plane.
        :param num_axis2_points: The number of y points to compute the fields for along the second axis in the plane.
        :param plane: One of {"xy", "xz", "yz"} to select which plane to take a cut of.
        :param offset: The offset from the edge of the geometry in the direction perpendicular to the plane to calc at.
        :return: Ex, Ey, and Ez, indexed by [z, y, x]. (One of these indices will only allow for the index 0).
        """
        print("Calculating field data")
        # Create a cuboid grid of points that the geometry is inscribed in
        x_points = [self.x_min+offset] if plane.upper() == "YZ" else np.linspace(self.x_min, self.x_max, num_axis1_points)
        y_points = [self.y_min+offset] if plane.upper() == "XZ" else np.linspace(self.y_min, self.y_max, num_axis1_points if plane.upper() == "YZ" else num_axis2_points)
        z_points = [self.z_min+offset] if plane.upper() == "XY" else np.linspace(self.z_min, self.z_max, num_axis2_points)
        num_x_points, num_y_points, num_z_points = len(x_points), len(y_points), len(z_points)
        Ex = np.zeros([num_z_points, num_y_points, num_x_points], dtype=complex)
        Ey = np.zeros([num_z_points, num_y_points, num_x_points], dtype=complex)
        Ez = np.zeros([num_z_points, num_y_points, num_x_points], dtype=complex)

        field_points = np.zeros([num_x_points * num_y_points * num_z_points, 3])
        # Iterate over the points
        for i in range(num_z_points):
            pt_z = z_points[i]
            for j in range(num_y_points):
                pt_y = y_points[j]
                for k in range(num_x_points):
                    pt_x = x_points[k]
                    field_points[k + j*num_x_points + i*num_x_points*num_y_points] = np.array([pt_x, pt_y, pt_z])

        # Find which tetrahedron each point lies in. This is vectorized.
        tet_indices = where(self.all_nodes, self.tets_node_ids, field_points)

        # Compute the field at each of the points
        for i, tet_index in enumerate(tet_indices):
            z_i = floor(i / (num_x_points * num_y_points)) % num_z_points
            y_i = floor(i / num_x_points) % num_y_points
            x_i = i % num_x_points
            # If a tetrahedron was not identified for the point, assign the field there to 0
            if tet_index == -1:
                Ex[z_i, y_i, x_i], Ey[z_i, y_i, x_i], Ez[z_i, y_i, x_i] = 0, 0, 0
                continue
            tet = self.tetrahedrons[tet_index]
            phis = [self.edge_coefficients[self.remap_edge_nums[edge]] if edge in self.remap_edge_nums else 0 for edge in tet.edges]
            # Note the indexing here is done with z_i first, y_i second, and x_i third. If we consider a 2D grid being
            # indexed, the first index corresponds to the row (vertical control), hence y_i second and x_i third.
            # Same idea applies to having z_i first.
            Ex[z_i, y_i, x_i], Ey[z_i, y_i, x_i], Ez[z_i, y_i, x_i] = tet.interpolate(phis, field_points[i])

        print("Finished calculating field data")
        return Ex, Ey, Ez

    def plot_fields(self, num_axis1_points=100, num_axis2_points=100, plane="xy", offset=0.1, phase=0., vmin=None, vmax=None, use_cached_fields=False, normalize=False):
        """
        Plot the fields in the selected plane. Note that field plotting is expensive due to needing to locate which
        tetrahedron each point lies in. Finer meshes may need to use fewer sample points.
        :param num_axis1_points: The number of points to compute the fields for along the first axis in the plane.
        :param num_axis2_points: The number of y points to compute the fields for along the second axis in the plane.
        :param plane: One of {"xy", "xz", "yz"} to select which plane to take a cut of.
        :param offset: The offset from the edge of the geometry in the direction perpendicular to the plane to calc at.
        :param phase: The phase in radians to calculate the fields at.
        :param vmin: The minimum value of the
        :param use_cached_fields: If ``True``, do not recompute the fields, just apply the phase and plot the result.
        :param normalize: If ``True``, normalize the maximum field values to 1. Note this causes a loss in the meaning
        of the relative magnitudes of Ex, Ey, Ez, since they are each separately normalized.
        :return: The figure containing all the field data (the result of running plt.figure()).
        """
        x_points = [self.x_min+offset] if plane.upper() == "YZ" else np.linspace(self.x_min, self.x_max, num_axis1_points)
        y_points = [self.y_min+offset] if plane.upper() == "XZ" else np.linspace(self.y_min, self.y_max, num_axis1_points if plane.upper() == "YZ" else num_axis2_points)
        z_points = [self.z_min+offset] if plane.upper() == "XY" else np.linspace(self.z_min, self.z_max, num_axis2_points)
        if not use_cached_fields or self.Ex is None or self.Ey is None or self.Ez is None:
            self.Ex, self.Ey, self.Ez = self.get_fields_in_plane(num_axis1_points, num_axis2_points, plane, offset)
        phase_shift = e**(-1j*phase)
        Ex, Ey, Ez = np.real(phase_shift*self.Ex), np.real(phase_shift*self.Ey), np.real(phase_shift*self.Ez)
        if normalize:
            Ex /= np.amax(np.abs(Ex))
            Ey /= np.amax(np.abs(Ey))
            Ez /= np.amax(np.abs(Ez))
        fig = plt.figure()
        plt.title(f"Fields in {plane.upper()}-plane, offset = {round(offset, 3)}")
        if plane.upper() == "YZ":
            axis1, axis2 = np.meshgrid(y_points, z_points)
            skip = (slice(None, None, 5), slice(None, None, 5))
            field_skip = (slice(None, None, 5), slice(None, None, 5), 0)
            plt.imshow(Ex[:, :, 0], extent=[self.y_min, self.y_max, self.z_min, self.z_max], cmap="cividis", vmin=vmin, vmax=vmax)
            plt.colorbar(label="Ex")
            plt.quiver(axis1[skip], axis2[skip], Ey[field_skip], Ez[field_skip], color="black")
        elif plane.upper() == "XZ":
            axis1, axis2 = np.meshgrid(x_points, z_points)
            skip = (slice(None, None, 5), slice(None, None, 5))
            field_skip = (slice(None, None, 5), 0, slice(None, None, 5))
            plt.imshow(Ey[:, 0, :], extent=[self.x_min, self.x_max, self.z_min, self.z_max], cmap="cividis", vmin=vmin, vmax=vmax)
            plt.colorbar(label="Ey")
            # plt.quiver(axis1[skip], axis2[skip], Ex[field_skip], Ez[field_skip], color="black")
        elif plane.upper() == "XY":
            axis1, axis2 = np.meshgrid(x_points, y_points)
            skip = (slice(None, None, 5), slice(None, None, 5))
            field_skip = (0, slice(None, None, 5), slice(None, None, 5))
            # Correct the orientation. TODO: Determine if this is necessary for the other planes
            Ez = np.flip(Ez, axis=1)
            plt.imshow(Ez[0, :, :], extent=[self.x_min, self.x_max, self.y_min, self.y_max], cmap="cividis", vmin=vmin, vmax=vmax)
            plt.colorbar(label="Ez")
            plt.quiver(axis1[skip], axis2[skip], Ex[field_skip], Ey[field_skip], color="black")
        else:
            raise RuntimeError(f"Invalid argument for plane '{plane}'. Must be one of 'xy', 'xz', or 'yz'.")
        return fig

    def plot_one_field(self, field, num_axis1_points=100, num_axis2_points=100, plane="xy", offset=0.1, phase=0., vmin=None, vmax=None, use_cached_fields=False, normalize=False):
        """
        Plot the fields in the selected plane. Note that field plotting is expensive due to needing to locate which
        tetrahedron each point lies in. Finer meshes may need to use fewer sample points.
        :param field: A string selecting which field to plot. One of "Ex", "Ey", or "Ez".
        :param num_axis1_points: The number of points to compute the fields for along the first axis in the plane.
        :param num_axis2_points: The number of y points to compute the fields for along the second axis in the plane.
        :param plane: One of {"xy", "xz", "yz"} to select which plane to take a cut of.
        :param offset: The offset from the edge of the geometry in the direction perpendicular to the plane to calc at.
        :param phase: The phase in radians to calculate the fields at.
        :param use_cached_fields: If ``True``, do not recompute the fields, just apply the phase and plot the result.
        :param normalize: If ``True``, normalize the maximum field values to 1. Note this causes a loss in the meaning
        of the relative magnitudes of Ex, Ey, Ez, since they are each separately normalized.
        :return: The figure containing all the field data (the result of running plt.figure()).
        """
        if not use_cached_fields or self.Ex is None or self.Ey is None or self.Ez is None:
            self.Ex, self.Ey, self.Ez = self.get_fields_in_plane(num_axis1_points, num_axis2_points, plane, offset)
        phase_shift = e**(-1j*phase)
        Ex, Ey, Ez = np.real(phase_shift*self.Ex), np.real(phase_shift*self.Ey), np.real(phase_shift*self.Ez)
        if field == "Ex":
            E_plot = Ex
        elif field == "Ey":
            E_plot = Ey
        elif field == "Ez":
            E_plot = Ez
        else:
            raise ValueError(f"Field {field} must be one of 'Ex', 'Ey', or 'Ez'.")
        if normalize:
            min, max = np.amin(E_plot), np.amax(E_plot)
            if abs(min) > abs(max):
                E_plot /= min
            else:
                E_plot /= max
            # E_plot /= np.amax(np.abs(E_plot))
        fig = plt.figure()
        plt.title(f"Fields in {plane.upper()}-plane, offset = {round(offset, 3)}")
        if plane.upper() == "YZ":
            plt.imshow(E_plot[:, :, 0], extent=[self.y_min, self.y_max, self.z_min, self.z_max], cmap="cividis", vmin=vmin, vmax=vmax)
            plt.colorbar(label=field)
        elif plane.upper() == "XZ":
            plt.imshow(E_plot[:, 0, :], extent=[self.x_min, self.x_max, self.z_min, self.z_max], cmap="cividis", vmin=vmin, vmax=vmax)
            plt.colorbar(label=field)
        elif plane.upper() == "XY":
            # Correct the orientation. TODO: Determine if this is necessary for the other planes
            E_plot = np.flip(E_plot, axis=1)
            plt.imshow(E_plot[0, :, :], extent=[self.x_min, self.x_max, self.y_min, self.y_max], cmap="cividis", vmin=vmin, vmax=vmax)
            plt.colorbar(label=field)
        else:
            raise RuntimeError(f"Invalid argument for plane '{plane}'. Must be one of 'xy', 'xz', or 'yz'.")
        return fig


count = 0


def gmres_cb(pr_norm):
    global count
    print(f"Iteration {count}")
    count += 1
    print(f"pr_norm = {pr_norm}")



