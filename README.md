# Navigation

This repository contains an inertial navigation system (INS) implementation for processing navigation data, attitude, and IMU measurements in both Python and MATLAB. It includes core coordinate-transform utilities, mechanization routines, Kalman filtering components, plotting helpers, and example scripts that demonstrate end-to-end use on real sensor datasets.

## Overview

The project is designed to estimate navigation states from inertial measurement data, including:

- position (latitude, longitude, height)
- velocity (North, East, Down)
- attitude (roll, pitch, yaw / quaternion)
- trajectory comparison against reference data

The implementation includes equivalent algorithms in both:

- Python under `src/python/`
- MATLAB under `src/MATLAB/`

## Repository structure

```text
Navigation/
├── data/                          # Navigation datasets used by examples
├── example/                      # Example scripts for Python and MATLAB
│   ├── main.py
│   ├── main_EKF.py
│   ├── main.m
│   └── main_EKF_dvl.m
├── src/
│   ├── python/
│   │   ├── __init__.py
│   │   ├── KalmanFilter.py
│   │   ├── SkewSymmetry.py
│   │   ├── calM_N.py
│   │   ├── coord_transform.py
│   │   ├── mechanization.py
│   │   ├── plotting.py
│   │   └── plotting_ekf.py
│   └── MATLAB/
│       ├── CoordTransform.m
│       ├── KalmanFilter.m
│       ├── Mechanization.m
│       ├── plotting.m
│       └── plotting_ekf.m
├── requirements.txt             # Python dependencies
├── LICENSE                      # if present in the repo, this is the project license
└── README.md                    # Project documentation
```

## Features

- INS mechanization using quaternion-based attitude propagation
- Coordinate transformation between Euler angles and quaternions
- Skew-symmetric matrix utilities
- Earth curvature / radii calculations (`calM_N`)
- Kalman filter implementations for navigation state estimation
- Visualization utilities for comparing reference and calculated trajectories
- Python and MATLAB parity for algorithm development and comparison

## Requirements

Python dependencies are listed in `requirements.txt` and include scientific computing libraries such as:

- NumPy
- SciPy
- Matplotlib

Install them with:

```bash
pip install -r requirements.txt
```

## Quick start (Python)

From the repository root, run the example script:

```bash
python example/main.py
```

This script:

1. loads the sample dataset from `data/path2sensor_output.mat`
2. runs the INS mechanization routine
3. compares estimated and reference states
4. plots the resulting navigation outputs

## Example usage

The Python implementation exposes an object-oriented mechanization workflow similar to:

```python
import os
import sys

THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.dirname(THIS_DIR)
SRC_DIR = os.path.join(ROOT_DIR, "src/python")

if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from mechanization import Mechanization_OOP

r_cal, v_cal, q_cal, roll_cal, pitch_cal, yaw_cal = Mechanization_OOP.runBatch(
    t, r_n, v_n, roll, pitch, yaw, f_b_IMU, w_b_ib_IMU
)
```

## Core modules

### `src/python/mechanization.py`
Contains the main INS mechanization logic, including single-step and batch propagation of position, velocity, and quaternion states.

### `src/python/coord_transform.py`
Provides quaternion and Euler conversion utilities used throughout the navigation pipeline.

### `src/python/KalmanFilter.py`
Implements filtering routines for state estimation and sensor fusion tasks.

### `src/python/plotting.py`
Generates plots for comparing estimated and reference navigation outputs.

### `src/MATLAB/`
Equivalent implementations in MATLAB for users working in the MATLAB ecosystem.

## Data

The repository includes navigation sensor data in the `data/` folder, which is used by the example scripts to validate the algorithms.

## Notes

This project is focused on navigation and inertial estimation workflows and is especially useful for:

- INS mechanization research
- attitude and trajectory estimation
- algorithm prototyping in Python or MATLAB
- comparing navigation results across frameworks

## Contributing

Contributions, bug fixes, and improvements are welcome. If you are extending the algorithm or adding new examples, please keep the Python and MATLAB implementations aligned where possible.

## License

Please check the repository for the project license file before using or distributing the code.

## Author

This repository appears to be authored by Arvin Chamasemani and is intended for navigation algorithm development and experimentation.
