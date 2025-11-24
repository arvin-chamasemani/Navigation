"""
INS Mechanization + EKF – Example Script (Python OOP Version)

This script demonstrates a full INS + EKF pipeline:
1) Load the MATLAB dataset
2) Run the OOP mechanization model
3) Run the error-state Kalman filter
4) Apply corrections to position / velocity / attitude
5) Plot reference vs. INS vs. EKF

Run from the examples/ folder:
    python run_ins_ekf_from_mat.py

programmer: Arvin Chamasemani
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
SRC_DIR = os.path.join(ROOT_DIR, "src", "python")

if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from mechanization import Mechanization_OOP         # provides run_batch(...)
from KalmanFilter import KalmanFilter          # provides run_batch_kf(...)
from coord_transform import CoordTransform      # provides q2mat, mat2euler
from plotting_ekf import plotting_ekf           # our INS vs EKF plotting
from SkewSymmetry import SkewSymmetry
from calM_N import calM_N


def main():
    # ---------------------------------------------------------
    # Load dataset
    # ---------------------------------------------------------
    data_path = os.path.join(THIS_DIR, "..", "data", "path2sensor_output.mat")
    data = loadmat(data_path)

    # Time and reference navigation states
    t    = data["t"].ravel()        # (N,)
    r_n  = data["r_n"]              # (3, N)
    v_n  = data["v_n"]              # (3, N)

    roll  = data["roll"].ravel()    # (N,)
    pitch = data["pitch"].ravel()   # (N,)
    yaw   = data["yaw"].ravel()     # (N,)

    # IMU inputs (noisy version)
    f_b_IMU    = data["f_b_IMU"]    # (3, N)
    w_b_ib_IMU = data["w_b_ib_IMU"] # (3, N)

    # Clean version (uncomment if you want ideal IMU)
    # f_b_IMU    = data["f_b"]
    # w_b_ib_IMU = data["w_b_ib"]

    # Sensor / measurement noise
    Q_sens  = data["Q"]             # (6, 6) sensor noise covariance
    var_GPS = data["var_GPS"].ravel()  # (3,) position meas. variance

    # ---------------------------------------------------------
    # Run INS mechanization
    # ---------------------------------------------------------
    print("Running INS mechanization...")

    # Expected output: r_cal, v_cal: (3, N), q_cal: (4, N), roll/pitch/yaw_cal: (N,)
    r_cal, v_cal, q_cal, roll_cal, pitch_cal, yaw_cal = Mechanization_OOP.runBatch(
        t, r_n, v_n, roll, pitch, yaw, f_b_IMU, w_b_ib_IMU
    )

    print("Mechanization completed.")

    # ---------------------------------------------------------
    # Run Kalman Filter
    # ---------------------------------------------------------
    print("Running error-state Kalman filter...")

    # Initial covariance P (9x9) for [dr; dv; dtheta]
    P0 = np.diag([1, 1, 1, 1e-3, 1e-3, 1e-3, 1e-9, 1e-9, 1e-9])

    x_hist, P_hist = KalmanFilter.run_batch_kf(
        t,
        r_cal, v_cal, q_cal,    # mechanization outputs
        r_n,                    # "measurements" (e.g. GPS position)
        f_b_IMU,                # IMU data
        Q_sens, var_GPS,        # noise covariances
        P0                      # initial covariance
    )

    print("Kalman filter completed.")

    # ---------------------------------------------------------
    # Apply EKF corrections to states
    # ---------------------------------------------------------
    print("Applying EKF corrections...")

    N = t.size
    r_ekf = np.zeros_like(r_cal)    # (3, N)
    v_ekf = np.zeros_like(v_cal)    # (3, N)

    roll_ekf  = np.zeros(N)
    pitch_ekf = np.zeros(N)
    yaw_ekf   = np.zeros(N)

    for k in range(N):
        # State error [dr; dv; dtheta]
        dr      = x_hist[0:3, k]
        dv      = x_hist[3:6, k]
        dtheta  = x_hist[6:9, k]

        # Position & velocity corrections
        r_ekf[:, k] = r_cal[:, k] + dr
        v_ekf[:, k] = v_cal[:, k] + dv

        # Attitude correction: C_bn_corr = (I + [dtheta]x) * C_bn
        C_bn_nom  = CoordTransform.q2mat(q_cal[:, k])
        C_bn_corr = (np.eye(3) + SkewSymmetry(dtheta)) @ C_bn_nom

        # Convert corrected DCM -> Euler (yaw, pitch, roll)
        yaw_k, pitch_k, roll_k = CoordTransform.mat2euler(C_bn_corr)
        yaw_ekf[k]   = yaw_k
        pitch_ekf[k] = pitch_k
        roll_ekf[k]  = roll_k

    print("EKF corrections applied.")

    # ---------------------------------------------------------
    # Plot reference vs INS vs EKF
    # ---------------------------------------------------------
    print("Plotting results (INS vs EKF)...")

    plotting_ekf(
        t,
        # reference
        r_n, v_n, roll, pitch, yaw,
        # INS (mechanization) outputs
        r_cal, v_cal, roll_cal, pitch_cal, yaw_cal,
        # EKF-corrected outputs
        r_ekf, v_ekf, roll_ekf, pitch_ekf, yaw_ekf,
    )

    print("Done.")


if __name__ == "__main__":
    main()
