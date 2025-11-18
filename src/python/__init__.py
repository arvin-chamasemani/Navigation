"""
Package initialization for the INS mechanization module.

This exposes the main public API:
- Mechanization_OOP   : OOP mechanization class
- CoordTransform      : Quaternion & Euler utilities
- calM_N              : Earth curvature radii
- SkewSymmetry        : Skew-symmetric matrix helper
- plot_INS            : Visualization utilities
"""

from .mechanization import Mechanization_OOP
from .coord_transform import CoordTransform
from .plotting import plot_INS

__all__ = [
    "Mechanization_OOP",
    "CoordTransform",
    "calM_N",
    "SkewSymmetry",
    "plot_INS",
]
