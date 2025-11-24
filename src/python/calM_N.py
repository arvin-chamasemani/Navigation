import numpy as np

def calM_N(phi):
    """
    Compute meridian (M) and prime-vertical (N) radii of curvature
    for a given geodetic latitude phi [rad] (IERS 2003 model).
    """
    a = 6378136.6  # semi-major axis [m]
    b = 6356751.9  # semi-minor axis [m]

    e = np.sqrt(1.0 - (b / a) ** 2)

    sin_phi = np.sin(phi)
    denom = np.sqrt(1.0 - (e * sin_phi) ** 2)

    N = a / denom
    M = a * (1.0 - e**2) / (denom**3)

    return M, N

