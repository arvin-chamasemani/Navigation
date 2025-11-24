import numpy as np
def SkewSymmetry(v):
    """
    Return the 3×3 skew-symmetric matrix S such that:
        S @ a = np.cross(v, a)
    for any 3×1 vector a.
    """
    v = np.asarray(v, dtype=float).reshape(3)
    return np.array(
        [
            [0.0, -v[2], v[1]],
            [v[2], 0.0, -v[0]],
            [-v[1], v[0], 0.0],
        ],
        dtype=float,
    )