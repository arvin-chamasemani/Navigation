"""
INS Mechanization – Example Script (Python OOP Version)

This script demonstrates a full INS processing pipeline:
1) Load the MATLAB dataset
2) Run the OOP mechanization model
3) Plot reference vs. estimated states

Run from the examples/ folder:
    python run_ins_from_mat.py
"""

import os
import sys
import numpy as np
from scipy.io import loadmat

# -------------------------------------------------------------
# Add src/ directory to Python path so project modules can be imported
# -------------------------------------------------------------
THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.dirname(THIS_DIR)
SRC_DIR = os.path.join(ROOT_DIR, "src/python")

if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from mechanization import Mechanization_OOP
from plotting import plot_INS

# -------------------------------------------------------------
# Load dataset
# -------------------------------------------------------------
data_path = os.path.join(THIS_DIR, "..", "data", "path2sensor_output.mat")
data = loadmat(data_path)

t    = data["t"].ravel()
r_n  = data["r_n"]
v_n  = data["v_n"]

roll  = data["roll"].ravel()
pitch = data["pitch"].ravel()
yaw   = data["yaw"].ravel()

# IMU inputs (choose noisy or clean)
f_b_IMU    = data["f_b_IMU"]
w_b_ib_IMU = data["w_b_ib_IMU"]

# Clean version:
# f_b_IMU    = data["f_b"]
# w_b_ib_IMU = data["w_b_ib"]

# -------------------------------------------------------------
# Run INS mechanization
# -------------------------------------------------------------
print("Running INS mechanization...")

r_cal, v_cal, q_cal, roll_cal, pitch_cal, yaw_cal = Mechanization_OOP.runBatch(
    t, r_n, v_n, roll, pitch, yaw, f_b_IMU, w_b_ib_IMU
)

print("Mechanization completed.")

# -------------------------------------------------------------
# Plot results
# -------------------------------------------------------------
print("Plotting results...")
plot_INS(
    t,
    r_n, v_n,
    roll, pitch, yaw,
    r_cal, v_cal,
    roll_cal, pitch_cal, yaw_cal,
)

print("Done.")
