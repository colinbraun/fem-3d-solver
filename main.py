from util import load_mesh, Edge, quad_eval, quad_sample_points, where
import numpy as np
from iwaveguide.waveguide import Waveguide
from scipy.linalg import inv
import matplotlib.pyplot as plt
from math import floor, e, pi

# Turn on interactive plotting
plt.ion()


class Waveguide3D:
    """Class representing a 3D waveguide with an input port and output port (not necessarily same size or shape)"""
    def __init__(self, filename):
        """
        Constructor of Waveguide3D.
        :param filename: A string giving the path to the .inp file to be loaded.
        """
        # Load the mesh
        self.all_nodes, self.tetrahedrons, self.tets_node_ids, self.all_edges, self.boundary_pec_edge_numbers, self.boundary_input_edge_numbers, self.boundary_output_edge_numbers, self.remap_edge_nums, self.all_edges_map, self.boundary_input_tets, self.boundary_output_tets = load_mesh(filename)
        # Find the minimum and maximum x, y, and z values of the waveguide.
        self.x_min = np.amin(self.all_nodes[:, 0])
        self.x_max = np.amax(self.all_nodes[:, 0])
        self.y_min = np.amin(self.all_nodes[:, 1])
        self.y_max = np.amax(self.all_nodes[:, 1])
        self.z_min = np.amin(self.all_nodes[:, 2])
        self.z_max = np.amax(self.all_nodes[:, 2])
        # Construct empty K and b matrices
        self.K = np.zeros([len(self.remap_edge_nums), len(self.remap_edge_nums)], dtype=complex)
        self.b = np.zeros([len(self.remap_edge_nums)], dtype=complex)
        # Create a Waveguide object of the input port
        self.input_port = Waveguide("rect_mesh_two_epsilons_coarse.inp", ["EB1", "EB2"], "EB3")
        # Set its mode to be the specified propgating mode (4 -> TM11)
        self.input_port.set_mode_index(4)
        # Solve the waveguide for k0 = 4
        self.input_port.solve_k0(4)
        # TODO: Change this to have an output port that differs from the input port. Will have to solve it too.
        self.output_port = self.input_port
        # Create an empty coefficient array. This is loaded by solve()
        self.edge_coefficients = np.zeros([len(self.remap_edge_nums), len(self.remap_edge_nums)])

    def solve(self):
        """
        Solve for the fields in the Waveguide3D object.
        :return: The edge coefficients to be applied to the interpolating functions (indexed by global edge number).
        """
        print("Begin constructing equation matrix")
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
                # Get the nodes that make up this edge
                node_il, node_jl = self.all_nodes[edge1.node1], self.all_nodes[edge1.node2]

                indices_l = [np.argwhere(tet.nodes == edge1.node1)[0][0], np.argwhere(tet.nodes == edge1.node2)[0][0]]
                # The simplex constants for nodes i and j of edge l
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
                # If we have an input port edge, we might need to perform the N_i integral (if this tet has a face on the port)
                # if edgei in boundary_input_edge_numbers:
                # If we are working with a tetrahedron with a face on the input port, need to perform integral with its edge
                if tet in self.boundary_input_tets:
                    # Need the third point that makes up the face of the input port
                    found_edge_nos = []
                    # Search for a different edge on the input port from this tetrahedral element
                    for edge in tet.edges:
                        if edge in self.boundary_input_edge_numbers:
                            # We found an edge containing the third node, note it
                            found_edge_nos.append(edge)
                        if len(found_edge_nos) == 2:
                            break
                    # It is possible we do not find another edge that lies on the input port. We are checking for this here.
                    if len(found_edge_nos) == 2:
                        found_edge1 = self.all_edges[found_edge_nos[0]]
                        found_edge2 = self.all_edges[found_edge_nos[1]]
                        nodes = np.unique(np.array([self.all_nodes[found_edge1.node1], self.all_nodes[found_edge1.node2], self.all_nodes[found_edge2.node1], self.all_nodes[found_edge2.node2]]), axis=0)
                        if len(nodes) != 3:
                            raise RuntimeError("Did not find 3 nodes on surface triangle")
                        # Perform the N_i \dot K_N integral for the surface using gaussian quadrature
                        sample_points = quad_sample_points(3, nodes[0], nodes[1], nodes[2])
                        N_i = np.zeros([len(sample_points), 3])
                        N_i[:, 0] = Axl + Bxl*sample_points[:, 1] + Cxl*sample_points[:, 2]
                        N_i[:, 1] = Ayl + Byl*sample_points[:, 0] + Cyl*sample_points[:, 2]
                        N_i[:, 2] = Azl + Bzl*sample_points[:, 0] + Czl*sample_points[:, 1]
                        # N_i = np.array([tet.interpolate([1] * 6, sample_point) for sample_point in sample_points])
                        # Compute the x component of the edge interpolating function for each of the sample points
                        # Get the E_inc field at each of the sample points
                        E_inc = np.array([np.array(self.input_port.get_field_at(sample_point[0], sample_point[1])) for sample_point in sample_points])
                        # Compute the dot product at each point
                        values = np.reshape(N_i[:, 0] * E_inc[:, 0] + N_i[:, 1] * E_inc[:, 1] + N_i[:, 2] * E_inc[:, 2], [len(sample_points), 1])
                        integral = quad_eval(nodes[0], nodes[1], nodes[2], values)
                        self.b[self.remap_edge_nums[self.all_edges_map[edge1]]] += -integral
                    else:
                        raise RuntimeError("Did not find 2 surface edges for tetrahedron marked as having a face on edge")

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
                        # TODO: Implement this integral and logic
                        # Perform the (n_hat x N_i) \dot (n_hat x N_j) integral
                        # For n_hat = z_hat, this turns out to be the same as N_i \dot N_j, so integral is same as inhomo wg
                        # Need all 3 nodes first
                        # Necessary constants from NASA paper eqs. 68-72.
                        # Need the third point that makes up the face of the input port
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
                            # TODO: Use the above properly (not being used at all right now for some reason)
                            ni_dot_nj = np.reshape(N_i[:, 0] * N_j[:, 0] + N_i[:, 1] * N_j[:, 1] + N_i[:, 2] * N_j[:, 2], [len(sample_points), 1])
                            # ni_dot_nj = N_i[:, 0] * N_j[:, 0] + N_i[:, 1] * N_j[:, 1] + N_i[:, 2] * N_j[:, 2]
                            integral = quad_eval(nodes[0], nodes[1], nodes[2], ni_dot_nj)
                            self.K[self.remap_edge_nums[self.all_edges_map[edge1]], self.remap_edge_nums[self.all_edges_map[edge2]]] += integral * self.input_port.get_selected_beta() * 1j

                    # If working with output port edges, need to do different integral
                    elif edgei in self.boundary_output_edge_numbers and edgej in self.boundary_output_edge_numbers:
                        # Perform the (n_hat x N_i) \dot (n_hat x N_j) integral
                        # For n_hat = z_hat, this turns out to be the same as N_i \dot N_j, so integral is same as inhomo wg
                        # Need all 3 nodes first
                        # Necessary constants from NASA paper eqs. 68-72.
                        # Need the third point that makes up the face of the input port
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
                            # TODO: Use the above properly (not being used at all right now for some reason)
                            ni_dot_nj = np.reshape(N_i[:, 0] * N_j[:, 0] + N_i[:, 1] * N_j[:, 1] + N_i[:, 2] * N_j[:, 2], [len(sample_points), 1])
                            integral = quad_eval(nodes[0], nodes[1], nodes[2], ni_dot_nj)
                            self.K[self.remap_edge_nums[self.all_edges_map[edge1]], self.remap_edge_nums[self.all_edges_map[edge2]]] += self.output_port.get_selected_beta() * integral * 1j
                    # Otherwise we are working with an inner edge (not on boundary). Do necessary integral.
                    # Currently removing this else unless it is found that it is needed
                    # else:
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
                    self.K[self.remap_edge_nums[self.all_edges_map[edge1]], self.remap_edge_nums[self.all_edges_map[edge2]]] += curl_dot_curl_part + dot_part
        print("Finished constructing equation matrix")
        print("Solving equation matrix")
        self.edge_coefficients = np.dot(inv(self.K), self.b)
        print("Finished solving equation matrix")
        return self.edge_coefficients

    def plot_fields(self, num_axis1_points=100, num_axis2_points=100, plane="xy", offset=0.1, phase=0.):
        """
        Plot the fields in the selected plane. Note that field plotting is expensive due to needing to locate which
        tetrahedron each point lies in. Finer meshes may need to use fewer sample points.
        :param num_axis1_points: The number of points to compute the fields for along the first axis in the plane.
        :param num_axis2_points: The number of y points to compute the fields for along the second axis in the plane.
        :param plane: One of {"xy", "xz", "yz"} to select which plane to take a cut of.
        :param offset: The offset from the edge of the geometry in the direction perpendicular to the plane to calc at.
        :param phase: The phase in radians to calculate the fields at.
        :return: The figure containing all the field data (the result of running plt.figure()).
        """
        print("Calculating field data")
        # Compute the bounds of the waveguide
        # Create a cuboid grid of points that the geometry is inscribed in
        x_points = [self.x_min+offset] if plane.upper() == "YZ" else np.linspace(self.x_min, self.x_max, num_axis1_points)
        y_points = [self.y_min+offset] if plane.upper() == "XZ" else np.linspace(self.y_min, self.y_max, num_axis1_points if plane.upper() == "YZ" else num_axis2_points)
        z_points = [self.z_min+offset] if plane.upper() == "XY" else np.linspace(self.z_min, self.z_max, num_axis2_points)
        num_x_points, num_y_points, num_z_points = len(x_points), len(y_points), len(z_points)
        Ex = np.zeros([num_z_points, num_y_points, num_x_points])
        Ey = np.zeros([num_z_points, num_y_points, num_x_points])
        Ez = np.zeros([num_z_points, num_y_points, num_x_points])

        field_points = np.zeros([num_x_points * num_y_points * num_z_points, 3])
        # Iterate over the points
        for i in range(num_z_points):
            pt_z = z_points[i]
            for j in range(num_y_points):
                pt_y = y_points[j]
                for k in range(num_x_points):
                    pt_x = x_points[k]
                    field_points[k + j*num_x_points + i*num_x_points*num_y_points] = np.array([pt_x, pt_y, pt_z])

        # Find which tetrahedron each point lies in
        tet_indices = where(self.all_nodes, self.tets_node_ids, field_points)

        shift = e**(1j*phase)
        # Compute the field at each of the points
        for i, tet_index in enumerate(tet_indices):
            tet = self.tetrahedrons[tet_index]
            phis = [self.edge_coefficients[self.remap_edge_nums[edge]] if edge in self.remap_edge_nums else 0 for edge in tet.edges]
            z_i = floor(i / (num_x_points * num_y_points)) % num_z_points
            y_i = floor(i / num_x_points) % num_y_points
            x_i = i % num_x_points
            # Note the indexing here is done with z_i first, y_i second, and x_i third. If we consider a 2D grid being
            # indexed, the first index corresponds to the row (vertical control), hence y_i second and x_i third.
            # Same idea applies to having z_i first.
            Ex[z_i, y_i, x_i], Ey[z_i, y_i, x_i], Ez[z_i, y_i, x_i] = np.real(np.multiply(tet.interpolate(phis, field_points[i]), shift))

        print("Finished calculating field data")

        fig = plt.figure()
        plt.title(f"Fields in {plane.upper()}-plane, offset = {round(offset, 3)}")
        if plane.upper() == "YZ":
            axis1, axis2 = np.meshgrid(y_points, z_points)
            skip = (slice(None, None, 5), slice(None, None, 5))
            field_skip = (slice(None, None, 5), slice(None, None, 5), 0)
            plt.imshow(Ex[:, :, 0], extent=[self.y_min, self.y_max, self.z_min, self.z_max], cmap="cividis")
            plt.colorbar(label="Ex")
            plt.quiver(axis1[skip], axis2[skip], Ey[field_skip], Ez[field_skip], color="black")
        elif plane.upper() == "XZ":
            axis1, axis2 = np.meshgrid(x_points, z_points)
            skip = (slice(None, None, 5), slice(None, None, 5))
            field_skip = (slice(None, None, 5), 0, slice(None, None, 5))
            plt.imshow(Ey[:, 0, :], extent=[self.x_min, self.x_max, self.z_min, self.z_max], cmap="cividis")
            plt.colorbar(label="Ey")
            plt.quiver(axis1[skip], axis2[skip], Ex[field_skip], Ez[field_skip], color="black")
        elif plane.upper() == "XY":
            axis1, axis2 = np.meshgrid(x_points, y_points)
            skip = (slice(None, None, 5), slice(None, None, 5))
            field_skip = (0, slice(None, None, 5), slice(None, None, 5))
            plt.imshow(Ez[0, :, :], extent=[self.x_min, self.x_max, self.y_min, self.y_max], cmap="cividis")
            plt.colorbar(label="Ez")
            plt.quiver(axis1[skip], axis2[skip], Ex[field_skip], Ey[field_skip], color="black")
        else:
            raise RuntimeError(f"Invalid argument for plane '{plane}'. Must be one of 'xy', 'xz', or 'yz'.")
        return fig


waveguide = Waveguide3D("rectangular_waveguide_3d_less_coarse.inp")
waveguide.solve()
for i in range(6):
    waveguide.plot_fields(plane="xy", offset=0.1, phase=i*pi/3)
