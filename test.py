# File to test out little pieces of code before using them in main.py
import numpy as np


def foo(x, y, z):
    """
    Test function for quadrature
    :param x: The x value to evaluate the function at.
    :param y: The y value to evaluate the function at.
    :param z: The z value to evaluate the function at.
    :return: The value of the function at a particular point.
    """
    return 3*x - 2*y


def quad_sample_points(n, p1, p2, p3):
    """
    Generate gaussian quadrature sample points for a triangle with points p1, p2, and p3.
    :param n: The order of the quadrature.
    :param p1: First point of a triangle.
    :param p2: Second point of a triangle.
    :param p3: Third point of a triangle.
    :return: An nx3 numpy array. Each row contains the x, y, and z values of a single sample point.
    """

    if n > 9 or n < 1:
        print("N must be a value between 1 and 9")
        return None

    loc = None
    x1, y1, z1 = p1[:]
    x2, y2, z2 = p2[:]
    x3, y3, z3 = p3[:]
    if n == 1:
        alpha = 1 / 3
        beta = 1 / 3
        temp1 = x1 * (1 - alpha - beta) + x2 * alpha + x3 * beta
        temp2 = y1 * (1 - alpha - beta) + y2 * alpha + y3 * beta
        temp3 = z1 * (1 - alpha - beta) + z2 * alpha + z3 * beta
        loc = [temp1, temp2, temp3]

    elif n == 2:
        alpha = [1 / 6, 2 / 3, 1 / 6]
        beta = [1 / 6, 1 / 6, 2 / 3]
        loc = np.zeros([3, 3])
        for m in range(len(alpha)):
            temp1 = x1 * (1 - alpha[m] - beta[m]) + x2 * alpha[m] + x3 * beta[m]
            temp2 = y1 * (1 - alpha[m] - beta[m]) + y2 * alpha[m] + y3 * beta[m]
            temp3 = z1 * (1 - alpha[m] - beta[m]) + z2 * alpha[m] + z3 * beta[m]
            loc[m, 0] = temp1
            loc[m, 1] = temp2
            loc[m, 2] = temp3
    elif n == 3:
        alpha = [1 / 3, 1 / 5, 1 / 5, 3 / 5]
        beta = [1 / 3, 3 / 5, 1 / 5, 1 / 5]
        loc = np.zeros([4, 3])
        for m in range(len(alpha)):
            temp1 = x1 * (1 - alpha[m] - beta[m]) + x2 * alpha[m] + x3 * beta[m]
            temp2 = y1 * (1 - alpha[m] - beta[m]) + y2 * alpha[m] + y3 * beta[m]
            temp3 = z1 * (1 - alpha[m] - beta[m]) + z2 * alpha[m] + z3 * beta[m]
            loc[m, 0] = temp1
            loc[m, 1] = temp2
            loc[m, 2] = temp3
    if loc is None:
        print("Error in computing quadrature. Did you choose a proper n value?")
    return loc


def quad_eval(p1, p2, p3, f):
    """
    Compute the integral using gaussian quadrature. The values of function (f(x, y, z)) should be stored in the argument f,
    stored as a column vector with each row corresponding to f(x, y) for a particular sample point. The sample points
    should come from quad_sample_points().
    :param p1: First point of a triangle.
    :param p2: Second point of a triangle.
    :param p3: Third point of a triangle.
    :param f: A numpy array containing f(x, y, z) in each row for each of the corresponding sample points.
    :return: The approximation of the integration (or exact value in some cases).
    """
    num_pts = len(f[:, 0])
    x1, y1, z1 = p1[:]
    x2, y2, z2 = p2[:]
    x3, y3, z3 = p3[:]

    if num_pts == 1:
        weights = [1]
    elif num_pts == 3:
        weights = [1/3,1/3,1/3]
    elif num_pts == 4:
        weights = [-27/48,25/48,25/48,25/48]
    elif num_pts == 6:
        weights = [0.22338158967801,0.22338158967801,0.22338158967801,
                   0.10995174365532,0.10995174365532,0.10995174365532]
    elif num_pts == 7:
        weights = [0.22500000000000,0.13239415278851,0.13239415278851,
                   0.13239415278851,0.12593918054483,0.12593918054483,
                   0.12593918054483]
    elif num_pts == 12:
        weights = [0.11678627572638,0.11678627572638,0.11678627572638,
                   0.05084490637021,0.05084490637021,0.05084490637021,
                   0.08285107561837,0.08285107561837,0.08285107561837,
                   0.08285107561837,0.08285107561837,0.08285107561837]
    elif num_pts == 13:
        weights = [-0.14957004446768,0.17561525743321,0.17561525743321,
                   0.17561525743321,0.05334723560884,0.05334723560884,
                   0.05334723560884,0.07711376089026,0.07711376089026,
                   0.07711376089026,0.07711376089026,0.07711376089026,
                   0.07711376089026]
    elif num_pts == 16:
        weights = [0.14431560767779,0.09509163426728,0.09509163426728,
                   0.09509163426728,0.10321737053472,0.10321737053472,
                   0.10321737053472,0.03245849762320,0.03245849762320,
                   0.03245849762320,0.02723031417443,0.02723031417443,
                   0.02723031417443,0.02723031417443,0.02723031417443,
                   0.02723031417443]
    elif num_pts == 19:
        weights = [0.097135796282799,0.031334700227139,0.031334700227139,
                   0.031334700227139,0.077827541004774,0.077827541004774,
                   0.077827541004774,0.079647738927210,0.079647738927210,
                   0.079647738927210,0.025577675658698,0.025577675658698,
                   0.025577675658698,0.043283539377289,0.043283539377289,
                   0.043283539377289,0.043283539377289,0.043283539377289,
                   0.043283539377289]
    else:
        print('Incorrect number of points used')

    val = 0
    Ak = np.sqrt(((y1-y3)*(z2-z3)-(y2-y3)*(z1-z3))**2 +
          ((x1-x3)*(z2-z3)-(x2-x3)*(z1-z3))**2 +
          ((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))**2)/2
    for m in range(len(weights)):
        val = val+weights[m]*f[m, :]
    val = val*Ak
    return val

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
p1 = [-1, 0, 0]
p2 = [1, 0, 0]
p3 = [0, 1, 0]
sample_points = quad_sample_points(n, p1, p2, p3)
values = np.zeros([len(sample_points), 1])
for i in range(len(sample_points)):
    values[i, 0] = foo(sample_points[i, 0], sample_points[i, 1], sample_points[i, 2])
result = quad_eval(p1, p2, p3, values)
