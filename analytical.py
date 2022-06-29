from math import sqrt, cos, sin
from scipy.constants import mu_0, epsilon_0, pi


def rect_wg_field_at(x, y, a, b, omega, m=1, n=0, permittivity=1., te=True):
    """
    Get the field at the specified point for a waveguide of size ``a`` x ``b`` with permittivity ``permittivity``.
    Assumes the bounds of the waveguide are [0, a] along x and [0, b] along y.
    :param x: The x coordinate to evaluate the field at.
    :param y: The y coordinate to evaluate the field at.
    :param a: The x length of the waveguide.
    :param b: The y length of the waveguide.
    :param omega: The angular frequency.
    :param m: The first mode index. Default = 1 (for TE10 mode)
    :param n: The second mode index. Default = 0 (for TE10 mode)
    :param permittivity: The permittivity of the material inside the waveguide. Default = 1
    :param te: If True, evaluates the field assuming a TE mode.
    :return: Ex, Ey, and Ez for the particular mode at the specified point.
    """
    k = omega * sqrt(mu_0*epsilon_0*permittivity)
    kc = sqrt((m*pi/a)**2 + (n*pi/b)**2)
    if k < kc:
        raise ValueError(f"Wave with k={k} is below cutoff frequency kc={kc}.")
    ex = omega * mu_0 * n * pi / kc**2 / b * cos(m*pi*x/a) * sin(n*pi*y/b)
    ey = -omega * mu_0 * m * pi / kc**2 / a * sin(m*pi*x/a) * cos(n*pi*y/b)
    return ex, ey, 0
