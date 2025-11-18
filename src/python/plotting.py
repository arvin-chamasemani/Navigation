import numpy as np
import matplotlib.pyplot as plt


def angle_error(a, b):
    """
    Wrapped angular difference: err = a - b in [-pi, pi]
    a, b: arrays (same shape) in radians
    """
    a = np.asarray(a)
    b = np.asarray(b)
    return np.arctan2(np.sin(a - b), np.cos(a - b))


def auto_error_ylim(y):
    """
    Simple automatic y-limits for error plots (Python version of AUTOERRORYLIM).
    """
    y = np.asarray(y).ravel()
    y = y[np.isfinite(y)]  # remove NaN / Inf

    if y.size == 0:
        return (-1.0, 1.0)

    n = y.size
    if n < 2:
        m = max(1.0, abs(y[0]))
        m *= 1.1
        return (-m, m)

    # Use middle 90% of data (ignore first/last 5%)
    i1 = max(0, int(np.floor(0.05 * n)))
    i2 = max(i1 + 1, int(np.ceil(0.95 * n)))
    s = y[i1:i2]

    s = s[np.isfinite(s)]
    if s.size == 0:
        s = y

    lo = np.min(s)
    hi = np.max(s)

    if lo == hi:
        span = max(1e-6, 0.1 * max(1.0, abs(lo)))
        return (lo - span, hi + span)
    else:
        pad = 0.1 * (hi - lo)
        return (lo - pad, hi + pad)


def plot_INS(t, r_n, v_n, roll, pitch, yaw,
             r_n_cal, v_n_cal, roll_cal, pitch_cal, yaw_cal):
    """
    Python version of the MATLAB plot_INS function.
    Compare reference vs. calculated INS states (velocity, position, attitude)
    and plot their errors.

    Parameters
    ----------
    t         : (N,)
    r_n       : (3, N) reference position [lat(rad), lon(rad), h(m)]
    v_n       : (3, N) reference velocity [N,E,D] (m/s)
    roll      : (N,)  reference roll (rad)
    pitch     : (N,)  reference pitch (rad)
    yaw       : (N,)  reference yaw (rad)
    r_n_cal   : (3, N) calculated position
    v_n_cal   : (3, N) calculated velocity
    roll_cal  : (N,)   calculated roll (rad)
    pitch_cal : (N,)   calculated pitch (rad)
    yaw_cal   : (N,)   calculated yaw (rad)
    """
    # Ensure arrays
    t = np.asarray(t).ravel()
    r_n = np.asarray(r_n)
    v_n = np.asarray(v_n)
    r_n_cal = np.asarray(r_n_cal)
    v_n_cal = np.asarray(v_n_cal)
    roll = np.asarray(roll).ravel()
    pitch = np.asarray(pitch).ravel()
    yaw = np.asarray(yaw).ravel()
    roll_cal = np.asarray(roll_cal).ravel()
    pitch_cal = np.asarray(pitch_cal).ravel()
    yaw_cal = np.asarray(yaw_cal).ravel()

    # ===================== VELOCITY (N,E,D) ============================
    vel_err = v_n_cal - v_n

    plt.figure(num='Velocity (NED) - Comparison & Error')
    # --- North velocity ---
    ax = plt.subplot(3, 2, 1)
    ax.plot(t, v_n_cal[0, :], linewidth=1.2)
    ax.plot(t, v_n[0, :], '--r', linewidth=1.0)
    ax.grid(True)
    ax.legend(['Calc', 'Ref'], loc='best')
    ax.set_title('North Velocity')
    ax.set_ylabel('v_N [m/s]')
    ax.set_xlabel('Time [s]')

    ax = plt.subplot(3, 2, 2)
    ax.plot(t, vel_err[0, :], 'k', linewidth=1.2)
    ax.grid(True)
    ax.set_title('North Velocity Error')
    ax.set_ylabel('Error [m/s]')
    ax.set_xlabel('Time [s]')
    ax.set_ylim(auto_error_ylim(vel_err[0, :]))

    # --- East velocity ---
    ax = plt.subplot(3, 2, 3)
    ax.plot(t, v_n_cal[1, :], linewidth=1.2)
    ax.plot(t, v_n[1, :], '--r', linewidth=1.0)
    ax.grid(True)
    ax.legend(['Calc', 'Ref'], loc='best')
    ax.set_title('East Velocity')
    ax.set_ylabel('v_E [m/s]')
    ax.set_xlabel('Time [s]')

    ax = plt.subplot(3, 2, 4)
    ax.plot(t, vel_err[1, :], 'k', linewidth=1.2)
    ax.grid(True)
    ax.set_title('East Velocity Error')
    ax.set_ylabel('Error [m/s]')
    ax.set_xlabel('Time [s]')
    ax.set_ylim(auto_error_ylim(vel_err[1, :]))

    # --- Down velocity ---
    ax = plt.subplot(3, 2, 5)
    ax.plot(t, v_n_cal[2, :], linewidth=1.2)
    ax.plot(t, v_n[2, :], '--r', linewidth=1.0)
    ax.grid(True)
    ax.legend(['Calc', 'Ref'], loc='best')
    ax.set_title('Down Velocity')
    ax.set_ylabel('v_D [m/s]')
    ax.set_xlabel('Time [s]')

    ax = plt.subplot(3, 2, 6)
    ax.plot(t, vel_err[2, :], 'k', linewidth=1.2)
    ax.grid(True)
    ax.set_title('Down Velocity Error')
    ax.set_ylabel('Error [m/s]')
    ax.set_xlabel('Time [s]')
    ax.set_ylim(auto_error_ylim(vel_err[2, :]))

    # ===================== POSITION (Lat, Lon, H) ======================
    lat_ref = np.rad2deg(r_n[0, :])
    lon_ref = np.rad2deg(r_n[1, :])
    h_ref = r_n[2, :]

    lat_cal = np.rad2deg(r_n_cal[0, :])
    lon_cal = np.rad2deg(r_n_cal[1, :])
    h_cal = r_n_cal[2, :]

    lat_err = lat_cal - lat_ref
    lon_err = lon_cal - lon_ref
    h_err = h_cal - h_ref

    plt.figure(num='Position - Comparison & Error')

    # --- Latitude ---
    ax = plt.subplot(3, 2, 1)
    ax.plot(t, lat_cal, linewidth=1.2)
    ax.plot(t, lat_ref, '--r', linewidth=1.0)
    ax.grid(True)
    ax.legend(['Calc', 'Ref'], loc='best')
    ax.set_title('Latitude')
    ax.set_ylabel('Latitude [deg]')
    ax.set_xlabel('Time [s]')

    ax = plt.subplot(3, 2, 2)
    ax.plot(t, lat_err, 'k', linewidth=1.2)
    ax.grid(True)
    ax.set_title('Latitude Error')
    ax.set_ylabel('Error [deg]')
    ax.set_xlabel('Time [s]')
    ax.set_ylim(auto_error_ylim(lat_err))

    # --- Longitude ---
    ax = plt.subplot(3, 2, 3)
    ax.plot(t, lon_cal, linewidth=1.2)
    ax.plot(t, lon_ref, '--r', linewidth=1.0)
    ax.grid(True)
    ax.legend(['Calc', 'Ref'], loc='best')
    ax.set_title('Longitude')
    ax.set_ylabel('Longitude [deg]')
    ax.set_xlabel('Time [s]')

    ax = plt.subplot(3, 2, 4)
    ax.plot(t, lon_err, 'k', linewidth=1.2)
    ax.grid(True)
    ax.set_title('Longitude Error')
    ax.set_ylabel('Error [deg]')
    ax.set_xlabel('Time [s]')
    ax.set_ylim(auto_error_ylim(lon_err))

    # --- Height ---
    ax = plt.subplot(3, 2, 5)
    ax.plot(t, h_cal, linewidth=1.2)
    ax.plot(t, h_ref, '--r', linewidth=1.0)
    ax.grid(True)
    ax.legend(['Calc', 'Ref'], loc='best')
    ax.set_title('Height')
    ax.set_ylabel('Height [m]')
    ax.set_xlabel('Time [s]')

    ax = plt.subplot(3, 2, 6)
    ax.plot(t, h_err, 'k', linewidth=1.2)
    ax.grid(True)
    ax.set_title('Height Error')
    ax.set_ylabel('Error [m]')
    ax.set_xlabel('Time [s]')
    ax.set_ylim(auto_error_ylim(h_err))

    # ===================== EULER ANGLES (Roll, Pitch, Yaw) =============
    roll_err = angle_error(roll_cal, roll)
    pitch_err = angle_error(pitch_cal, pitch)
    yaw_err = angle_error(yaw_cal, yaw)

    roll_ref_deg = np.rad2deg(roll)
    pitch_ref_deg = np.rad2deg(pitch)
    yaw_ref_deg = np.rad2deg(yaw)

    roll_cal_deg = np.rad2deg(roll_cal)
    pitch_cal_deg = np.rad2deg(pitch_cal)
    yaw_cal_deg = np.rad2deg(yaw_cal)

    roll_err_deg = np.rad2deg(roll_err)
    pitch_err_deg = np.rad2deg(pitch_err)
    yaw_err_deg = np.rad2deg(yaw_err)

    plt.figure(num='Euler Angles - Comparison & Error')

    # --- Roll ---
    ax = plt.subplot(3, 2, 1)
    ax.plot(t, roll_cal_deg, linewidth=1.2)
    ax.plot(t, roll_ref_deg, '--r', linewidth=1.0)
    ax.grid(True)
    ax.legend(['Calc', 'Ref'], loc='best')
    ax.set_title('Roll')
    ax.set_ylabel('[deg]')
    ax.set_xlabel('Time [s]')

    ax = plt.subplot(3, 2, 2)
    ax.plot(t, roll_err_deg, 'k', linewidth=1.2)
    ax.grid(True)
    ax.set_title('Roll Error')
    ax.set_ylabel('Error [deg]')
    ax.set_xlabel('Time [s]')
    ax.set_ylim(auto_error_ylim(roll_err_deg))

    # --- Pitch ---
    ax = plt.subplot(3, 2, 3)
    ax.plot(t, pitch_cal_deg, linewidth=1.2)
    ax.plot(t, pitch_ref_deg, '--r', linewidth=1.0)
    ax.grid(True)
    ax.legend(['Calc', 'Ref'], loc='best')
    ax.set_title('Pitch')
    ax.set_ylabel('[deg]')
    ax.set_xlabel('Time [s]')

    ax = plt.subplot(3, 2, 4)
    ax.plot(t, pitch_err_deg, 'k', linewidth=1.2)
    ax.grid(True)
    ax.set_title('Pitch Error')
    ax.set_ylabel('Error [deg]')
    ax.set_xlabel('Time [s]')
    ax.set_ylim(auto_error_ylim(pitch_err_deg))

    # --- Yaw ---
    ax = plt.subplot(3, 2, 5)
    ax.plot(t, yaw_cal_deg, linewidth=1.2)
    ax.plot(t, yaw_ref_deg, '--r', linewidth=1.0)
    ax.grid(True)
    ax.legend(['Calc', 'Ref'], loc='best')
    ax.set_title('Yaw')
    ax.set_ylabel('[deg]')
    ax.set_xlabel('Time [s]')

    ax = plt.subplot(3, 2, 6)
    ax.plot(t, yaw_err_deg, 'k', linewidth=1.2)
    ax.grid(True)
    ax.set_title('Yaw Error')
    ax.set_ylabel('Error [deg]')
    ax.set_xlabel('Time [s]')
    ax.set_ylim(auto_error_ylim(yaw_err_deg))

    plt.tight_layout()
    plt.show()
