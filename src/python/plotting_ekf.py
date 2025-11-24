import numpy as np
import matplotlib.pyplot as plt


def angle_error(a, b):
    """
    Wrapped angular difference: err = a - b in [-pi, pi]
    a, b: arrays of same shape (radians)
    """
    a = np.asarray(a)
    b = np.asarray(b)
    return np.arctan2(np.sin(a - b), np.cos(a - b))


def auto_error_ylim(y):
    """
    Simple automatic y-limits for error plots.
    Follows the same logic as the MATLAB version.
    """
    y = np.asarray(y).ravel()
    y = y[np.isfinite(y)]

    if y.size == 0:
        return (-1.0, 1.0)

    n = y.size
    if n < 2:
        m = max(1.0, abs(y[0]))
        span = 1.1 * m
        return (-span, span)

    # use middle 90% (ignore first/last 5%)
    i1 = max(0, int(np.floor(0.05 * n)))
    i2 = max(i1 + 1, int(np.ceil(0.95 * n)))
    s = y[i1:i2]
    s = s[np.isfinite(s)]
    if s.size == 0:
        s = y

    lo = s.min()
    hi = s.max()

    if lo == hi:
        span = max(1e-6, 0.1 * max(1.0, abs(lo)))
        return (lo - span, hi + span)
    else:
        pad = 0.1 * (hi - lo)
        return (lo - pad, hi + pad)


def plotting_ekf(
    t,
    r_n, v_n, roll, pitch, yaw,                # reference
    r_n_cal, v_n_cal, roll_cal, pitch_cal, yaw_cal,  # mechanization
    r_n_ekf, v_n_ekf, roll_ekf, pitch_ekf, yaw_ekf   # EKF
):
    """
    Compare reference vs mechanization vs EKF (velocity, position, attitude)
    and plot both mechanization and EKF errors on the same axes.
    All angles are assumed in radians.
    """

    t = np.asarray(t)

    # ===================== VELOCITY (N,E,D) ============================
    vel_err_mech = v_n_cal - v_n
    vel_err_ekf  = v_n_ekf - v_n

    plt.figure(figsize=(10, 8), num="Velocity (NED) - Comparison & Error")

    # --- North velocity ---
    plt.subplot(3, 2, 1)
    plt.plot(t, v_n_cal[0, :], "b",  linewidth=1.2, label="INS")
    plt.plot(t, v_n_ekf[0, :], "g--", linewidth=1.2, label="EKF")
    plt.plot(t, v_n[0, :],     "r:",  linewidth=1.2, label="Ref")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("North Velocity")
    plt.ylabel("v_N [m/s]")
    plt.xlabel("Time [s]")

    plt.subplot(3, 2, 2)
    plt.plot(t, vel_err_mech[0, :], "k", linewidth=1.2, label="INS Error")
    plt.plot(t, vel_err_ekf[0, :],  "m", linewidth=1.2, label="EKF Error")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("North Velocity Error")
    plt.ylabel("Error [m/s]")
    plt.xlabel("Time [s]")
    plt.ylim(auto_error_ylim(np.concatenate([vel_err_mech[0, :], vel_err_ekf[0, :]])))

    # --- East velocity ---
    plt.subplot(3, 2, 3)
    plt.plot(t, v_n_cal[1, :], "b",  linewidth=1.2, label="INS")
    plt.plot(t, v_n_ekf[1, :], "g--", linewidth=1.2, label="EKF")
    plt.plot(t, v_n[1, :],     "r:",  linewidth=1.2, label="Ref")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("East Velocity")
    plt.ylabel("v_E [m/s]")
    plt.xlabel("Time [s]")

    plt.subplot(3, 2, 4)
    plt.plot(t, vel_err_mech[1, :], "k", linewidth=1.2, label="INS Error")
    plt.plot(t, vel_err_ekf[1, :],  "m", linewidth=1.2, label="EKF Error")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("East Velocity Error")
    plt.ylabel("Error [m/s]")
    plt.xlabel("Time [s]")
    plt.ylim(auto_error_ylim(np.concatenate([vel_err_mech[1, :], vel_err_ekf[1, :]])))

    # --- Down velocity ---
    plt.subplot(3, 2, 5)
    plt.plot(t, v_n_cal[2, :], "b",  linewidth=1.2, label="INS")
    plt.plot(t, v_n_ekf[2, :], "g--", linewidth=1.2, label="EKF")
    plt.plot(t, v_n[2, :],     "r:",  linewidth=1.2, label="Ref")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("Down Velocity")
    plt.ylabel("v_D [m/s]")
    plt.xlabel("Time [s]")

    plt.subplot(3, 2, 6)
    plt.plot(t, vel_err_mech[2, :], "k", linewidth=1.2, label="INS Error")
    plt.plot(t, vel_err_ekf[2, :],  "m", linewidth=1.2, label="EKF Error")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("Down Velocity Error")
    plt.ylabel("Error [m/s]")
    plt.xlabel("Time [s]")
    plt.ylim(auto_error_ylim(np.concatenate([vel_err_mech[2, :], vel_err_ekf[2, :]])))

    # ===================== POSITION (Lat, Lon, H) ======================

    lat_ref = np.rad2deg(r_n[0, :])
    lon_ref = np.rad2deg(r_n[1, :])
    h_ref   = r_n[2, :]

    lat_mech = np.rad2deg(r_n_cal[0, :])
    lon_mech = np.rad2deg(r_n_cal[1, :])
    h_mech   = r_n_cal[2, :]

    lat_ekf = np.rad2deg(r_n_ekf[0, :])
    lon_ekf = np.rad2deg(r_n_ekf[1, :])
    h_ekf   = r_n_ekf[2, :]

    lat_err_mech = lat_mech - lat_ref
    lon_err_mech = lon_mech - lon_ref
    h_err_mech   = h_mech   - h_ref

    lat_err_ekf = lat_ekf - lat_ref
    lon_err_ekf = lon_ekf - lon_ref
    h_err_ekf   = h_ekf   - h_ref

    plt.figure(figsize=(10, 8), num="Position - Comparison & Error")

    # --- Latitude ---
    plt.subplot(3, 2, 1)
    plt.plot(t, lat_mech, "b",  linewidth=1.2, label="INS")
    plt.plot(t, lat_ekf,  "g--", linewidth=1.2, label="EKF")
    plt.plot(t, lat_ref,  "r:",  linewidth=1.2, label="Ref")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("Latitude")
    plt.ylabel("Latitude [deg]")
    plt.xlabel("Time [s]")

    plt.subplot(3, 2, 2)
    plt.plot(t, lat_err_mech, "k", linewidth=1.2, label="INS Error")
    plt.plot(t, lat_err_ekf,  "m", linewidth=1.2, label="EKF Error")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("Latitude Error")
    plt.ylabel("Error [deg]")
    plt.xlabel("Time [s]")
    plt.ylim(auto_error_ylim(np.concatenate([lat_err_mech, lat_err_ekf])))

    # --- Longitude ---
    plt.subplot(3, 2, 3)
    plt.plot(t, lon_mech, "b",  linewidth=1.2, label="INS")
    plt.plot(t, lon_ekf,  "g--", linewidth=1.2, label="EKF")
    plt.plot(t, lon_ref,  "r:",  linewidth=1.2, label="Ref")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("Longitude")
    plt.ylabel("Longitude [deg]")
    plt.xlabel("Time [s]")

    plt.subplot(3, 2, 4)
    plt.plot(t, lon_err_mech, "k", linewidth=1.2, label="INS Error")
    plt.plot(t, lon_err_ekf,  "m", linewidth=1.2, label="EKF Error")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("Longitude Error")
    plt.ylabel("Error [deg]")
    plt.xlabel("Time [s]")
    plt.ylim(auto_error_ylim(np.concatenate([lon_err_mech, lon_err_ekf])))

    # --- Height ---
    plt.subplot(3, 2, 5)
    plt.plot(t, h_mech, "b",  linewidth=1.2, label="INS")
    plt.plot(t, h_ekf,  "g--", linewidth=1.2, label="EKF")
    plt.plot(t, h_ref,  "r:",  linewidth=1.2, label="Ref")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("Height")
    plt.ylabel("Height [m]")
    plt.xlabel("Time [s]")

    plt.subplot(3, 2, 6)
    plt.plot(t, h_err_mech, "k", linewidth=1.2, label="INS Error")
    plt.plot(t, h_err_ekf,  "m", linewidth=1.2, label="EKF Error")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("Height Error")
    plt.ylabel("Error [m]")
    plt.xlabel("Time [s]")
    plt.ylim(auto_error_ylim(np.concatenate([h_err_mech, h_err_ekf])))

    # ===================== EULER ANGLES (Roll, Pitch, Yaw) =============

    roll_err_mech  = angle_error(roll_cal,  roll)
    pitch_err_mech = angle_error(pitch_cal, pitch)
    yaw_err_mech   = angle_error(yaw_cal,   yaw)

    roll_err_ekf  = angle_error(roll_ekf,  roll)
    pitch_err_ekf = angle_error(pitch_ekf, pitch)
    yaw_err_ekf   = angle_error(yaw_ekf,   yaw)

    roll_ref_deg  = np.rad2deg(roll)
    pitch_ref_deg = np.rad2deg(pitch)
    yaw_ref_deg   = np.rad2deg(yaw)

    roll_mech_deg  = np.rad2deg(roll_cal)
    pitch_mech_deg = np.rad2deg(pitch_cal)
    yaw_mech_deg   = np.rad2deg(yaw_cal)

    roll_ekf_deg  = np.rad2deg(roll_ekf)
    pitch_ekf_deg = np.rad2deg(pitch_ekf)
    yaw_ekf_deg   = np.rad2deg(yaw_ekf)

    roll_err_mech_deg  = np.rad2deg(roll_err_mech)
    pitch_err_mech_deg = np.rad2deg(pitch_err_mech)
    yaw_err_mech_deg   = np.rad2deg(yaw_err_mech)

    roll_err_ekf_deg  = np.rad2deg(roll_err_ekf)
    pitch_err_ekf_deg = np.rad2deg(pitch_err_ekf)
    yaw_err_ekf_deg   = np.rad2deg(yaw_err_ekf)

    plt.figure(figsize=(10, 8), num="Euler Angles - Comparison & Error")

    # --- Roll ---
    plt.subplot(3, 2, 1)
    plt.plot(t, roll_mech_deg, "b",  linewidth=1.2, label="INS")
    plt.plot(t, roll_ekf_deg,  "g--", linewidth=1.2, label="EKF")
    plt.plot(t, roll_ref_deg,  "r:",  linewidth=1.2, label="Ref")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("Roll")
    plt.ylabel("[deg]")
    plt.xlabel("Time [s]")

    plt.subplot(3, 2, 2)
    plt.plot(t, roll_err_mech_deg, "k", linewidth=1.2, label="INS Error")
    plt.plot(t, roll_err_ekf_deg,  "m", linewidth=1.2, label="EKF Error")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("Roll Error")
    plt.ylabel("Error [deg]")
    plt.xlabel("Time [s]")
    plt.ylim(auto_error_ylim(np.concatenate([roll_err_mech_deg, roll_err_ekf_deg])))

    # --- Pitch ---
    plt.subplot(3, 2, 3)
    plt.plot(t, pitch_mech_deg, "b",  linewidth=1.2, label="INS")
    plt.plot(t, pitch_ekf_deg,  "g--", linewidth=1.2, label="EKF")
    plt.plot(t, pitch_ref_deg,  "r:",  linewidth=1.2, label="Ref")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("Pitch")
    plt.ylabel("[deg]")
    plt.xlabel("Time [s]")

    plt.subplot(3, 2, 4)
    plt.plot(t, pitch_err_mech_deg, "k", linewidth=1.2, label="INS Error")
    plt.plot(t, pitch_err_ekf_deg,  "m", linewidth=1.2, label="EKF Error")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("Pitch Error")
    plt.ylabel("Error [deg]")
    plt.xlabel("Time [s]")
    plt.ylim(auto_error_ylim(np.concatenate([pitch_err_mech_deg, pitch_err_ekf_deg])))

    # --- Yaw ---
    plt.subplot(3, 2, 5)
    plt.plot(t, yaw_mech_deg, "b",  linewidth=1.2, label="INS")
    plt.plot(t, yaw_ekf_deg,  "g--", linewidth=1.2, label="EKF")
    plt.plot(t, yaw_ref_deg,  "r:",  linewidth=1.2, label="Ref")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("Yaw")
    plt.ylabel("[deg]")
    plt.xlabel("Time [s]")

    plt.subplot(3, 2, 6)
    plt.plot(t, yaw_err_mech_deg, "k", linewidth=1.2, label="INS Error")
    plt.plot(t, yaw_err_ekf_deg,  "m", linewidth=1.2, label="EKF Error")
    plt.grid(True)
    plt.legend(loc="best")
    plt.title("Yaw Error")
    plt.ylabel("Error [deg]")
    plt.xlabel("Time [s]")
    plt.ylim(auto_error_ylim(np.concatenate([yaw_err_mech_deg, yaw_err_ekf_deg])))

    plt.tight_layout()
    plt.show()
