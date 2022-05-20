# File to test out little pieces of code before using them in main.py
import numpy as np
from iwaveguide.waveguide import Waveguide
from util import quad_eval, quad_sample_points


def foo(x, y, z):
    """
    Test function for quadrature
    :param x: The x value to evaluate the function at.
    :param y: The y value to evaluate the function at.
    :param z: The z value to evaluate the function at.
    :return: The value of the function at a particular point.
    """
    return 2 - x - 2*y


x1, y1, z1 = 1, 2, 3
x2, y2, z2 = 1, 2, 3
x3, y3, z3 = 1, 2, 3
x4, y4, z4 = 1, 2, 3
init_array = np.array([[1, x1, y1, z1],
                       [1, x2, y2, z2],
                       [1, x3, y3, z3],
                       [1, x4, y4, z4]], dtype=float)

all_cofactors = np.zeros([4, 4])
# Iterate over each row
for row in range(4):
    cofactors = np.zeros([4])
    # Iterate over each column, computing the cofactor determinant of the row + column combination
    for col in range(4):
        # Compute the cofactor (remove the proper row and column and compute the determinant)
        cofactors[col] = np.linalg.det(np.delete(np.delete(init_array, row, axis=0), col, axis=1))
    all_cofactors[row] = cofactors

n = 2
# p1 = [-1, 0, 0]
# p2 = [1, 0, 0]
# p3 = [0, 1, 0]
p1 = [0, 0, 0]
p2 = [1, 1/2, 0]
p3 = [0, 1, 0]
sample_points = quad_sample_points(n, p1, p2, p3)
values = np.zeros([len(sample_points), 1])
for i in range(len(sample_points)):
    values[i, 0] = foo(sample_points[i, 0], sample_points[i, 1], sample_points[i, 2])
# The result of this should be 1/3 (and has been verified to be this)
result = quad_eval(p1, p2, p3, values)

# waveguide_2d = Waveguide("rect_mesh_two_epsilons_coarse.inp", 2, [1, 1])

# betas, all_eigenvectors, k0s = waveguide_2d.solve()
# waveguide_2d.plot_dispersion(k0s, betas)
