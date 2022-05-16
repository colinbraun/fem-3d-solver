import numpy as np
import math


class Edge:
    """Class representing an edge. Immutable (Values should only be read from, not written to)."""
    def __init__(self, node1, node2):
        """
        Constructor
        :param node1: The global node number of the first node of the edge
        :param node2: The global node number of the second node of the edge
        """
        self.node1 = node1
        self.node2 = node2
        # print(len(Element.all_nodes))
        node1_t, node2_t = TriangleElement.all_nodes[node1], TriangleElement.all_nodes[node2]
        # TODO: Evaluate whether the sign should be considered or not (currently does assign it a sign)
        sign_multiplier = 1 if node1 < node2 else -1
        self.length = sign_multiplier * math.sqrt((node2_t[0] - node1_t[0])**2 + (node2_t[1] - node1_t[1])**2)

    def flip(self):
        """Return a copy of this object with the node numbers switched"""
        return Edge(self.node2, self.node1)

    def line(self):
        """Return a tuple of matrices, the first containing a list of 2 x points, the second of y points"""
        return [TriangleElement.all_nodes[self.node1][0], TriangleElement.all_nodes[self.node2][0]], [TriangleElement.all_nodes[self.node1][1], TriangleElement.all_nodes[self.node2][1]]

    def __eq__(self, other):
        """Edges are considered equal if they have the same global node numbers"""
        return self.node1 == other.node1 and self.node2 == other.node2 or self.node1 == other.node2 and self.node2 == other.node1

    def __hash__(self):
        """
        Hash the edge object for fast performance on maps/sets.
        Two edges with flipped node numbers produce the same hash value (are considered equal when compared thru map/set)
        """
        if self.node1 < self.node2:
            hash_value = (self.node1, self.node2).__hash__()
        else:
            hash_value = (self.node2, self.node1).__hash__()
        # print(f"Hash value: {hash_value}")
        return hash_value


class TriangleElement:
    """A class representing a triangular element (containing 3 global node numbers and 3 global edge numbers, as well as a permittivity) """
    all_nodes = []
    all_edges = []

    def __init__(self, nodes, edges, permittivity=1):
        """
        :param nodes: Three global node numbers as a numpy array
        :param edges: Three global edge numbers as a numpy array
        :param permittivity: The relative permittivity associated with this element (relative to vacuum permittivity)
        """
        self.nodes = nodes
        self.edges = edges
        self.permittivity = permittivity

    def area(self):
        """
        Compute the area of this element
        :return: The area of the triangle element
        """
        x1, y1 = TriangleElement.all_nodes[self.nodes[0]]
        x2, y2 = TriangleElement.all_nodes[self.nodes[1]]
        x3, y3 = TriangleElement.all_nodes[self.nodes[2]]
        return area(x1, y1, x2, y2, x3, y3)

    def is_inside(self, x, y):
        """
        Check if a point is inside this triangle element
        :param x: The x coordinate of the point
        :param y: The y coordinate of the point
        :return: True if the point lies in the triangle element, false otherwise
        """
        node1, node2, node3 = TriangleElement.all_nodes[self.nodes[0]], TriangleElement.all_nodes[self.nodes[1]], TriangleElement.all_nodes[self.nodes[2]]
        return is_inside(node1, node2, node3, x, y)

    def is_adjacent_to(self, element):
        """
        Check if this element is adjacent to the passed element
        :param element: The element to check if this is adjacent to
        :return: True if it is, false otherwise
        """
        node_count = 0
        if self.nodes[0] in element.nodes:
            node_count += 1
        if self.nodes[1] in element.nodes:
            node_count += 1
        if self.nodes[2] in element.nodes:
            node_count += 1
        return node_count == 2

    def nodal_interpolate(self, phi1, phi2, phi3, x, y):
        x1, y1 = TriangleElement.all_nodes[self.nodes[0]][0], TriangleElement.all_nodes[self.nodes[0]][1]
        x2, y2 = TriangleElement.all_nodes[self.nodes[1]][0], TriangleElement.all_nodes[self.nodes[1]][1]
        x3, y3 = TriangleElement.all_nodes[self.nodes[2]][0], TriangleElement.all_nodes[self.nodes[2]][1]
        a1 = x2 * y3 - x3 * y2
        a2 = x3 * y1 - x1 * y3
        a3 = x1 * y2 - x2 * y1
        b1 = y2 - y3
        b2 = y3 - y1
        b3 = y1 - y2
        c1 = x3 - x2
        c2 = x1 - x3
        c3 = x2 - x1

        n1 = phi1 * (a1 + b1 * x + c1 * y)
        n2 = phi2 * (a2 + b2 * x + c2 * y)
        n3 = phi3 * (a3 + b3 * x + c3 * y)
        return (n1 + n2 + n3) / 2 / self.area()

    def edge_interpolate(self, phi1, phi2, phi3, x, y):
        area = self.area()
        x_component = 0
        y_component = 0
        phis = [phi1, phi2, phi3]
        count = 0
        for edge_number in self.edges:
            # Generate the edge from 2 of the nodes (done in a CCW fashion by choice of tuples in for loop)
            edge = TriangleElement.all_edges[edge_number]
            # TODO: Fix the below statement so that it doesn't print TRUE
            # if edge != Element.all_edges[self.edges[count]]:
            #     print("TRUE")
            # Index of the third node (the one not making up the edge) going in a CCW fashion
            node1, node2 = edge.node1, edge.node2
            node3 = set(self.nodes).difference({node1, node2}).pop()
            n1_index = np.where(self.nodes == node1)[0][0]
            n2_index = np.where(self.nodes == node2)[0][0]
            # n2_index = self.nodes.index(node2)
            # print(str(n1_index[0][0]) + str(n2_index[0][0]))
            negate = 1
            match str(n1_index) + str(n2_index):
                case "01":
                    n3_index = 2
                case "02":
                    n1_index, n2_index, n3_index, negate = 2, 0, 1, -1
                case "10":
                    n1_index, n2_index, n3_index, negate = 0, 1, 2, -1
                case "12":
                    n1_index, n2_index, n3_index = 1, 2, 0
                case "20":
                    n1_index, n2_index, n3_index = 2, 0, 1
                case "21":
                    n1_index, n2_index, n3_index, negate = 1, 2, 0, -1

            # Create the ccw node list started from the first node of the edge
            nodes_lk = (TriangleElement.all_nodes[self.nodes[n1_index]], TriangleElement.all_nodes[self.nodes[n2_index]], TriangleElement.all_nodes[self.nodes[n3_index]])

            a_i_l, a_j_l = nodes_lk[1][0] * nodes_lk[2][1] - nodes_lk[2][0] * nodes_lk[1][1], nodes_lk[2][0] * nodes_lk[0][1] - nodes_lk[0][0] * nodes_lk[2][1]
            b_i_l, b_j_l = nodes_lk[1][1] - nodes_lk[2][1], nodes_lk[2][1] - nodes_lk[0][1]
            c_i_l, c_j_l = nodes_lk[2][0] - nodes_lk[1][0], nodes_lk[0][0] - nodes_lk[2][0]

            A_l = a_i_l * b_j_l - a_j_l * b_i_l
            B_l = c_i_l * b_j_l - c_j_l * b_i_l
            C_l = a_i_l * c_j_l - a_j_l * c_i_l
            D_l = b_i_l * c_j_l - b_j_l * c_i_l

            x_component += phis[count] * negate * edge.length / 4 / area**2 * (A_l + B_l*y)
            y_component += phis[count] * negate * edge.length / 4 / area**2 * (C_l + D_l*x)
            count += 1
        return x_component, y_component

    def __eq__(self, other):
        if self.nodes[0] in other.nodes and self.nodes[1] in other.nodes and self.nodes[2] in other.nodes:
            return True
        return False


class TetrahedralElement:
    """Class representing a tetrahedral (3D) element consisting of 4 TriangularElement objects"""

    def __init__(self, edges, permittivity=1):
        """

        :param edges: A numpy array containing the 6 edges (as global edge numbers) of the tetrahedral element
        :param permittivity: The permittivity of this element
        """
        self.edges = edges
        self.permittivity = permittivity
        # Store the unique set of points that make up this tetrahedron in no particular order
        self.points = TriangleElement.all_nodes[np.unique([edge.node1 for edge in [TriangleElement.all_edges[i] for i in edges]] + [edge.node2 for edge in [TriangleElement.all_edges[i] for i in edges]])]
        # print("Filler")

    def volume(self):
        """
        Compute and return the volume of this tetrahedral element
        :return: The volume of this tetrahedral element
        """
        # TODO: Consider just making this a field of the class
        # Volume calculated using eq (157) in NASA paper
        mat = [[1, self.points[0][0], self.points[0][1], self.points[0][2]],
               [1, self.points[1][0], self.points[1][1], self.points[1][2]],
               [1, self.points[2][0], self.points[2][1], self.points[2][2]],
               [1, self.points[3][0], self.points[3][1], self.points[3][2]]]
        return abs(np.linalg.det(mat) / 6)


def load_mesh_block(filename, block_name):
    """
    A function to load a block of data from a .inp file (abaqus format).
    Note that the node or element numbers are chopped off from the returned data.
    Node numbers contained in blocks of elements should be used as indices into the ALLNODES block.
    :param filename: The file to load the block of mesh data from
    :param block_name: The name of the block in the file that should be loaded
    :return: A 2D numpy array, each row containing a node or element (element can be 2D or 3D)
    """

    with open(filename, 'r') as file:
        lines = file.readlines()
        count = 0
        line = lines[count]
        # Go until we get to named section
        while block_name not in line:
            count += 1
            line = lines[count]
        # Skip over the line with the block_name in it
        count += 1
        line = lines[count]
        skip = count
        while "*" not in line:
            count += 1
            line = lines[count]
        rows = count - skip
        # Load all the node information into a numpy array, slicing off the node numbers (these are implied by index)
        # If we are dealing with elements, subtract off 1 so that the global node numbers work as indices
        if block_name != "ALLNODES":
            data = np.subtract(np.loadtxt(filename, skiprows=skip, max_rows=rows, delimiter=',', dtype=np.int32)[:, 1:], 1)
        else:
            data = np.loadtxt(filename, skiprows=skip, max_rows=rows, delimiter=',')[:, 1:]
        return data


def load_mesh(filename):
    """
    Load a mesh from a file. Must have at least 2 blocks, one containing all surface element_to_node_conn, the other the edge ones.
    :param filename: The name of the mesh file (a .inp a.k.a. abaqus file) to load
    :return: A number of items. See comments at bottom of this function, above the return statement.
    """
    # Load the nodes
    all_nodes = load_mesh_block(filename, "ALLNODES")
    # Make the nodes accessible globally (TODO: Do this differently?)
    TriangleElement.all_nodes = all_nodes
    # Load the tetrahedrons
    element_to_node_conn = load_mesh_block(filename, "Tetrahedrons")
    # All the edges (the index is the global edge number, i.e. all_edges[0] gets edge w/ global edge number 0)
    all_edges = []
    TriangleElement.all_edges = all_edges
    # A temporary map from an edge object to its global edge number (MUCH lower computational complexity)
    all_edges_map = {}
    # Keep track of what edge number we are on
    edge_count = 0

    # Create a list of tetrahedrons
    tetrahedrons = []
    # Iterate over all the tetrahedrons
    for tet in element_to_node_conn:
        # The list to hold the 4 triangle elements that make up the tetrahedron
        triangle_elements = []
        # The set holding the 6 global edge numbers of this tetrahedron
        tet_edges = set()
        # Iterate over each triangle face of the tetrahedron
        for element in (tet.take((0, 1, 2)), tet.take((0, 1, 3)), tet.take((0, 2, 3)), tet.take((1, 2, 3))):
            # Construct 3 edges
            edge1 = Edge(element[0], element[1])
            edge2 = Edge(element[0], element[2])
            edge3 = Edge(element[1], element[2])
            # For each of the 3 edges in this element, check if we have created this edge before. If not, create it.
            for edge in (edge1, edge2, edge3):
                if edge in all_edges_map:
                    # We do not want duplicate edges in our "all_edges" (global edge) list
                    pass
                else:
                    # Otherwise, we have not come across this edge before, so add it
                    all_edges.append(edge)
                    all_edges_map[edge] = edge_count
                    edge_count += 1
            # Get the global edge numbers for these 3 edges
            edge1_number = all_edges_map[edge1]
            edge2_number = all_edges_map[edge2]
            edge3_number = all_edges_map[edge3]
            tet_edges.add(edge1_number)
            tet_edges.add(edge2_number)
            tet_edges.add(edge3_number)
            # triangle = TriangleElement(np.array(element), np.array((edge1_number, edge2_number, edge3_number)), 1)
            # triangle_elements.append(triangle)
        # Create the tetrahedron object, converting the global edge numbers set to a numpy array
        tetrahedrons.append(TetrahedralElement(np.array(np.array(list(tet_edges)))))

    # Convert tetrahedrons to numpy array
    all_tets = np.array(tetrahedrons)
    # Convert edges to numpy array
    all_edges = np.array(all_edges)
    # Make it globally available
    TriangleElement.all_edges = all_edges

    # Load the PEC Wall triangle elements
    boundary_pec_elements = load_mesh_block(filename, "PECWalls")
    boundary_pec_edges = [Edge(element[0], element[1]) for element in boundary_pec_elements]
    boundary_pec_edge_numbers = set(all_edges_map[edge] for edge in boundary_pec_edges)
    # Load the InputPort triangle elements
    boundary_input_elements = load_mesh_block(filename, "InputPort")
    boundary_input_edges = [Edge(element[0], element[1]) for element in boundary_input_elements]
    boundary_input_edge_numbers = set(all_edges_map[edge] for edge in boundary_input_edges)
    # Load the OutputPort triangle elements
    boundary_output_elements = load_mesh_block(filename, "OutputPort")
    boundary_output_edges = [Edge(element[0], element[1]) for element in boundary_output_elements]
    boundary_output_edge_numbers = set(all_edges_map[edge] for edge in boundary_output_edges)

    # Get the set of non-boundary global edge numbers
    inner_edge_numbers = set(np.arange(0, len(all_edges))) - boundary_pec_edge_numbers - boundary_input_edge_numbers - boundary_output_edge_numbers
    # A map that takes one of the inner edge numbers and maps it to a unique integer between [0, number of inner edges]
    remap_inner_edge_nums = {item: i for i, item in enumerate(inner_edge_numbers)}

    # Return the results of loading the mesh
    # all_tets: A list of all the TetrahedralElement objects that make up the entire geometry
    # all_edges: A list of all the global edge numbers that make up the entire geometry
    # boundary_pec_edge_numbers: A set of all the global edge numbers that lie on the PEC wall of the geometry
    # boundary_input_edge_numbers: A set of all the global edge numbers that lie on the InputPort wall of the geometry
    # boundary_output_edge_numbers: A set of all the global edge numbers that lie on the OutputPort wall of the geometry
    # remap_inner_edge_nums: A map that takes one of the inner edge numbers and maps it to a unique integer between [0, number of inner edges]
    # all_edges_map: A map from an Edge object to its global edge number
    return all_tets, all_edges, boundary_pec_edge_numbers, boundary_input_edge_numbers, boundary_output_edge_numbers, remap_inner_edge_nums, all_edges_map


def area(x1, y1, x2, y2, x3, y3):
    return abs((x1 * (y2 - y3) + x2 * (y3 - y1)
                + x3 * (y1 - y2)) / 2.0)


def sign(p1, p2, p3):
    return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])


# A function to check whether point P(x, y)
# lies inside the triangle formed by
# A(x1, y1), B(x2, y2) and C(x3, y3)
def is_inside(pt1, pt2, pt3, x, y):
    if sign((x, y), pt1, pt2) < 0.0:
        return False
    if sign((x, y), pt2, pt3) < 0.0:
        return False
    if sign((x, y), pt3, pt1) < 0.0:
        return False
    return True
