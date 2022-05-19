from util import load_mesh_block
from util import load_mesh
import numpy as np

# print(load_mesh_block("rectangular_waveguide_3d.inp", "ALLNODES"))
# print(load_mesh_block("rectangular_waveguide_3d.inp", "InputPort"))
all_nodes, tetrahedrons, all_edges, boundary_pec_edge_numbers, boundary_input_edge_numbers, boundary_output_edge_numbers, remap_edge_nums, all_edges_map = load_mesh("rectangular_waveguide_3d.inp")
# Initialize the K and b matrices
K = np.zeros([len(remap_edge_nums), len(remap_edge_nums)])
b = np.zeros([len(remap_edge_nums)])

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

