from util import load_mesh, Edge
import numpy as np
from iwaveguide.waveguide import Waveguide

# print(load_mesh_block("rectangular_waveguide_3d.inp", "ALLNODES"))
# print(load_mesh_block("rectangular_waveguide_3d.inp", "InputPort"))
all_nodes, tetrahedrons, all_edges, boundary_pec_edge_numbers, boundary_input_edge_numbers, boundary_output_edge_numbers, remap_edge_nums, all_edges_map, boundary_input_triangles, boundary_output_triangles = load_mesh("rectangular_waveguide_3d.inp")
# Initialize the K and b matrices
K = np.zeros([len(remap_edge_nums), len(remap_edge_nums)], dtype=complex)
b = np.zeros([len(remap_edge_nums)], dtype=complex)

# TODO: Get the proper beta value from the inhomogeneous waveguide code
beta = 1

# Iterate over the tetrahedrons and construct the K and b matrices (this concept comes from Jin page 454)
for tet in tetrahedrons:

    # Compute the x-mean, y-mean, and z-mean of the points that make up the tetrahedron (needed later)
    x_mean = np.average(tet.points[:, 0])
    y_mean = np.average(tet.points[:, 1])
    z_mean = np.average(tet.points[:, 2])
    # Start by finding the simplex (barycentric) coordinates for the nodes
    # Each row is for a node. Each column is for a, b, c, and d (in order) from NASA paper eq. 162
    all_cofactors = np.zeros([4, 4])
    # Iterate over each row
    for row in range(4):
        cofactors = np.zeros([4])
        # Iterate over each column, computing the cofactor determinant of the row + column combination
        for col in range(4):
            # Compute the cofactor (remove the proper row and column and compute the determinant)
            cofactors[col] = np.linalg.det(np.delete(np.delete(np.append(tet.points, np.ones([4, 1]), 1), row, axis=0), col, axis=1))
        all_cofactors[row] = cofactors

    # Iterate over the edges of the tetrahedron
    for edgei in tet.edges:
        # Skip over PEC walls
        if edgei in boundary_pec_edge_numbers:
            continue
        # Get a hold of the Edge object
        edge1 = all_edges[edgei]
        # Get the nodes that make up this edge
        node_il, node_jl = all_nodes[edge1.node1], all_nodes[edge1.node2]

        indices_l = [np.argwhere(tet.nodes == edge1.node1)[0][0], np.argwhere(tet.nodes == edge1.node2)[0][0]]
        # The simplex coordinates for nodes i and j of edge l
        # Necessary constants from NASA paper eqs. 163-172
        a_il, a_jl = all_cofactors[indices_l]
        Axl = a_il[0]*a_jl[1] - a_il[1]*a_jl[0]
        Bxl = a_il[2]*a_jl[1] - a_il[1]*a_jl[2]
        Cxl = a_il[3]*a_jl[2] - a_il[2]*a_jl[3]
        Ayl = a_il[0]*a_jl[2] - a_il[2]*a_jl[0]
        Byl = a_il[1]*a_jl[2] - a_il[2]*a_jl[1]
        Cyl = a_il[3]*a_jl[2] - a_il[2]*a_jl[3]
        Azl = a_il[0]*a_jl[3] - a_il[3]*a_jl[0]
        Bzl = a_il[1]*a_jl[3] - a_il[3]*a_jl[1]
        Czl = a_il[2]*a_jl[3] - a_il[3]*a_jl[2]
        # Iterate over the edges of the tetrahedron
        for edgej in tet.edges:
            # Skip over PEC walls
            if edgej in boundary_pec_edge_numbers:
                continue
            edge2 = all_edges[edgej]
            node_ik, node_jk = all_nodes[edge1.node1], all_nodes[edge1.node2]

            # Find the indices of the edge of interest
            indices_k = [np.argwhere(tet.nodes == edge2.node1)[0][0], np.argwhere(tet.nodes == edge2.node2)[0][0]]
            # The simplex coordinates for nodes i and j of edge l
            a_ik, a_jk = all_cofactors[indices_k]
            # Necessary constants from NASA paper eqs. 163-172
            Axk = a_ik[0] * a_jk[1] - a_ik[1] * a_jk[0]
            Bxk = a_ik[2] * a_jk[1] - a_ik[1] * a_jk[2]
            Cxk = a_ik[3] * a_jk[2] - a_ik[2] * a_jk[3]
            Ayk = a_ik[0] * a_jk[2] - a_ik[2] * a_jk[0]
            Byk = a_ik[1] * a_jk[2] - a_ik[2] * a_jk[1]
            Cyk = a_ik[3] * a_jk[2] - a_ik[2] * a_jk[3]
            Azk = a_ik[0] * a_jk[3] - a_ik[3] * a_jk[0]
            Bzk = a_ik[1] * a_jk[3] - a_ik[3] * a_jk[1]
            Czk = a_ik[2] * a_jk[3] - a_ik[3] * a_jk[2]
            # If working with input port edges, need to do different integral
            if edgei in boundary_input_edge_numbers and edgej in boundary_input_edge_numbers:
                # TODO: Implement this integral and logic
                # Perform the (n_hat x N_i) \dot (n_hat x N_j) integral
                # For n_hat = z_hat, this turns out to be the same as N_i \dot N_j, so integral is same as inhomo wg
                # Need all 3 nodes first
                # Necessary constants from NASA paper eqs. 68-72.
                pass
            # If working with output port edges, need to do different integral
            elif edgei in boundary_output_edge_numbers and edgej in boundary_output_edge_numbers:
                # TODO: Implement this integral and logic
                pass
            # Otherwise we are working with an inner edge (not on boundary). Do necessary integral.
            else:
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
                dot_part = tet.permittivity * edge1.length * edge2.length / 1296 / tet.volume**3 * i_sum

# Get a Waveguide object
waveguide_2d = Waveguide("iwaveguide/rect_mesh_two_epsilons_coarse.inp")

# Iterate over the input port TriangleElement objects
for element in boundary_input_triangles:
    nodes = element.nodes
    area = element.area()
    # Iterate over the 3 edges of the element
    for l, k in ((0, 1), (1, 2), (2, 0)):
        # The first edge whose basis vector will be integrated against another edge (edge2)
        edge1 = Edge(nodes[l], nodes[k])
        # Skip edges on the boundary
        if all_edges_map[edge1] in boundary_input_edge_numbers:
            continue
        # Index of the third node (the one not making up the edge) going in a CCW fashion
        m = (k + 1) % 3
        # Create the ccw node list started from the first node of the edge
        nodes_lk = (all_nodes[nodes[l]], all_nodes[nodes[k]], all_nodes[nodes[m]])
        # These constants are for the first edge (call it edge l)
        a_i_l, a_j_l = nodes_lk[1][0] * nodes_lk[2][1] - nodes_lk[2][0] * nodes_lk[1][1], nodes_lk[2][0] * nodes_lk[0][
            1] - nodes_lk[0][0] * nodes_lk[2][1]
        b_i_l, b_j_l = nodes_lk[1][1] - nodes_lk[2][1], nodes_lk[2][1] - nodes_lk[0][1]
        c_i_l, c_j_l = nodes_lk[2][0] - nodes_lk[1][0], nodes_lk[0][0] - nodes_lk[2][0]
        # These come from the edge basis function definitions (see NASA paper eqs. 51 and 56-60)
        A_l = a_i_l * b_j_l - a_j_l * b_i_l
        B_l = c_i_l * b_j_l - c_j_l * b_i_l
        C_l = a_i_l * c_j_l - a_j_l * c_i_l
        D_l = b_i_l * c_j_l - b_j_l * c_i_l

        # Perform the N_i \dot K_N integral

        # Iterate over the 3 edges of the element (p stands for prime here, see Jin's book pg 455)
        for l_p, k_p in ((0, 1), (1, 2), (2, 0)):
            # The second edge whose basis vector will be integrated against the other edge (edge1)
            edge2 = Edge(nodes[l_p], nodes[k_p])
            # Skip edges on the boundary
            if all_edges_map[edge2] in boundary_input_edge_numbers:
                continue
            # Index of the third node (the one not making up the edge) going in a CCW fashion
            m_p = (k_p + 1) % 3
            # Create the ccw node list started from the first node of the edge
            nodes_lpkp = (all_nodes[nodes[l_p]], all_nodes[nodes[k_p]], all_nodes[nodes[m_p]])
            # Constants needed to calculate the integral involving the 2 edges
            # These come from the nodal basis function definitions (see Jin's book pg 441)
            # These constants are for the second edge (call it edge k)
            a_i_k, a_j_k = nodes_lpkp[1][0] * nodes_lpkp[2][1] - nodes_lpkp[2][0] * nodes_lpkp[1][1], nodes_lpkp[2][0] * \
                           nodes_lpkp[0][1] - nodes_lpkp[0][0] * nodes_lpkp[2][1]
            b_i_k, b_j_k = nodes_lpkp[1][1] - nodes_lpkp[2][1], nodes_lpkp[2][1] - nodes_lpkp[0][1]
            c_i_k, c_j_k = nodes_lpkp[2][0] - nodes_lpkp[1][0], nodes_lpkp[0][0] - nodes_lpkp[2][0]
            # These come from the edge basis function definitions (see NASA paper eqs. 51 and 56-60)
            A_k = a_i_k * b_j_k - a_j_k * b_i_k
            B_k = c_i_k * b_j_k - c_j_k * b_i_k
            C_k = a_i_k * c_j_k - a_j_k * c_i_k
            D_k = b_i_k * c_j_k - b_j_k * c_i_k
            # See An Introduction to the finite element method 3rd edition page 426 and NASA paper eqs. 73-77
            I1 = A_l * A_k + C_l * C_k
            I2 = (C_l * D_k + C_k * D_l) * x_mean
            I3 = (A_l * B_k + A_k * B_l) * y_mean
            I4 = B_l * B_k * (sum([node[1] ** 2 for node in nodes_lk]) + 9 * y_mean ** 2) / 12
            I5 = D_l * D_k * (sum([node[0] ** 2 for node in nodes_lk]) + 9 * x_mean ** 2) / 12
            # Ni_bar dot Nj_bar - NASA paper pg 12
            value = edge1.length * edge2.length * (I1 + I2 + I3 + I4 + I5) / 16 / (area ** 3)
            K[remap_edge_nums[all_edges_map[edge1]], remap_edge_nums[all_edges_map[edge2]]] += value * beta * 1j
