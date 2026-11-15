# Navigation

A navigation and inertial estimation project focused on modeling, filtering, and visualizing motion using inertial measurement data. This repository brings together Python and MATLAB implementations of a full inertial navigation system (INS), including coordinate transformations, mechanization, Kalman filtering, and trajectory plotting.

## About this project

This project was built to explore how inertial navigation algorithms can be implemented and validated from raw IMU-like data. It combines theory and practical computation in a compact, readable codebase that is useful for research, prototyping, and algorithm comparison across two popular engineering environments: Python and MATLAB.

The repository emphasizes:

- inertial navigation state estimation
- position and velocity propagation
- attitude estimation using quaternions
- coordinate transformations between Euler angles and quaternions
- comparison of estimated vs. reference trajectories
- clean visualization for analysis and presentation

## Project highlights

- Python implementation for modern scientific computing workflows
- MATLAB implementation for engineering and research pipelines
- Core INS mechanization logic in a reusable OOP structure
- Rotation and transformation utilities for navigation mathematics
- Plotting tools for analyzing position, velocity, and orientation
- Example scripts demonstrating full end-to-end processing

## Repository structure

```text
Navigation/
├── data/                     # Sample navigation datasets
├── example/                  # Example scripts for Python and MATLAB
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
├── requirements.txt
├── README.md
└── LICENSE
```

## Core capabilities

### INS mechanization
The project includes a mechanization pipeline that propagates:

- position in latitude, longitude, and altitude
- velocity in NED coordinates
- attitude represented by quaternions
- roll, pitch, and yaw estimates

### Coordinate transformations
The library includes quaternion and Euler angle conversion utilities, enabling robust manipulation of orientation states throughout the navigation workflow.

### Filtering and estimation
Kalman filter modules provide the groundwork for state estimation and filtering of inertial sensor data, supporting navigation estimation tasks.

### Visualization
The plotting utilities generate plots for comparing reference trajectories against estimated navigation states, making it easy to validate algorithm performance.

## Getting started

### Prerequisites

This project uses Python scientific packages. Install dependencies with:

```bash
pip install -r requirements.txt
```

### Run the Python example

```bash
python example/main.py
```

This example loads the dataset in `data/`, runs the mechanization pipeline, and produces plots comparing estimated and reference navigation states.

## Example usage

```python
from mechanization import Mechanization_OOP

r_cal, v_cal, q_cal, roll_cal, pitch_cal, yaw_cal = Mechanization_OOP.runBatch(
    t, r_n, v_n, roll, pitch, yaw, f_b_IMU, w_b_ib_IMU
)
```

## Technologies used

- Python
- NumPy
- SciPy
- Matplotlib
- MATLAB

## Why this project matters

Navigation systems are at the heart of autonomous vehicles, aerospace, robotics, and sensor fusion applications. This project showcases how foundational navigation equations and practical estimation methods can be implemented in a clear and reusable way, making it a strong example of applied computational engineering.

## Future directions

Possible next steps include:

- improved sensor fusion with GPS and IMU integration
- extended EKF/UKF-based filtering
- robust error modeling and bias estimation
- package structure improvements for easier importing
- documentation and benchmarking for algorithm comparison

## Author

Arvin Chamasemani

## License

This project is provided under the repository's license terms. Please check the repository for the relevant license file before reuse or redistribution.
