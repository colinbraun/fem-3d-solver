from util import load_mesh_block
from util import load_mesh
import numpy as np

# print(load_mesh_block("rectangular_waveguide_3d.inp", "ALLNODES"))
# print(load_mesh_block("rectangular_waveguide_3d.inp", "InputPort"))
tetrahedrons, all_edges, boundary_pec_edge_numbers, boundary_input_edge_numbers, boundary_output_edge_numbers, remap_inner_edge_nums, all_edges_map = load_mesh("rectangular_waveguide_3d.inp")
# print(tetrahedrons[0].faces[0].nodes[0])
# TODO: Replace this temporary placeholder with a real map from the mesh-loading code
Att = np.zeros([len(remap_inner_edge_nums), len(remap_inner_edge_nums)])
k = np.zeros([len(remap_inner_edge_nums)])

# Iterate over the tetrahedrons and construct the K and b matrices (this concept comes from Jin page 454)
for tet in tetrahedrons:
    # print(tet.edges)
    # print(tet.volume())
    # Iterate over the edges of the tetrahedron
    for edgei in tet.edges:
        # Iterate over the edges of the tetrahedron
        for edgej in tet.edges:
            pass

