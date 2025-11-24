import numpy as np
from coord_transform import CoordTransform
from SkewSymmetry import SkewSymmetry
from calM_N import calM_N
# programmer: Arvin chamasrmani

class Mechanization_OOP:
    """
    OOP-style INS mechanization.

    Main usage:
        r_cal, v_cal, q_cal, roll_cal, pitch_cal, yaw_cal = Mechanization_OOP.runBatch(...)
    """

    @staticmethod
    def runBatch(t, r_n, v_n, roll, pitch, yaw, f_b, w_b_ib):
        """
        Run the INS mechanization over a full time sequence.

        Parameters
        ----------
        t       : (N,)    time vector [s]
        r_n     : (3, N)  reference position [lat, lon, h]
        v_n     : (3, N)  reference velocity [N, E, D] (m/s)
        roll    : (N,)    initial roll sequence (rad)
        pitch   : (N,)    initial pitch sequence (rad)
        yaw     : (N,)    initial yaw sequence (rad)
        f_b     : (3, N)  specific force in body frame (m/s^2)
        w_b_ib  : (3, N)  body angular rate (rad/s)

        Returns
        -------
        r_n_cal : (3, N)  estimated position
        v_n_cal : (3, N)  estimated velocity
        q_cal   : (4, N)  estimated quaternion states
        roll_cal, pitch_cal, yaw_cal : (N,) Euler angles from q_cal
        """
        t = np.asarray(t)
        r_n = np.asarray(r_n)
        v_n = np.asarray(v_n)
        roll = np.asarray(roll)
        pitch = np.asarray(pitch)
        yaw = np.asarray(yaw)
        f_b = np.asarray(f_b)
        w_b_ib = np.asarray(w_b_ib)

        N = t.size

        # Initial states from first sample
        r_n_0 = r_n[:, 0]
        v_n_0 = v_n[:, 0]
        q_0 = CoordTransform.euler2q(roll[0], pitch[0], yaw[0])

        r_n_cal = np.zeros_like(r_n)
        v_n_cal = np.zeros_like(v_n)
        q_cal = np.zeros((4, N), dtype=float)

        r_n_cal[:, 0] = r_n_0
        v_n_cal[:, 0] = v_n_0
        q_cal[:, 0] = q_0

        # Time stepping (keeps same fixed dt as original MATLAB code)
        for i in range(1, N):
            # Using a constant step: dt = t(5) - t(4) in MATLAB
            dt_i = t[4] - t[3]

            r_n_cal[:, i], v_n_cal[:, i], q_cal[:, i] = Mechanization_OOP.dynamics(
                r_n_cal[:, i - 1],
                v_n_cal[:, i - 1],
                q_cal[:, i - 1],
                f_b[:, i - 1],
                w_b_ib[:, i - 1],
                dt_i,
            )

        # Convert quaternions back to Euler angles for analysis
        yaw_cal = np.zeros(N, dtype=float)
        pitch_cal = np.zeros(N, dtype=float)
        roll_cal = np.zeros(N, dtype=float)

        for i in range(N):
            yaw_i, pitch_i, roll_i = CoordTransform.q2euler(q_cal[:, i])
            yaw_cal[i] = yaw_i
            pitch_cal[i] = pitch_i
            roll_cal[i] = roll_i

        return r_n_cal, v_n_cal, q_cal, roll_cal, pitch_cal, yaw_cal

    @staticmethod
    def dynamics(r_n, v_n, q, f_b, omega_ib_b, dt):
        """
        Single-step INS update.

        Performs:
            1) attitude update (quaternion)
            2) velocity update (in NED)
            3) position update (lat, lon, h)

        Parameters
        ----------
        r_n        : (3,) position [lat, lon, h]
        v_n        : (3,) velocity [N, E, D]
        q          : (4,) attitude quaternion [x, y, z, w]
        f_b        : (3,) specific force in body frame
        omega_ib_b : (3,) body angular rate
        dt         : float, time step [s]

        Returns
        -------
        r_n_out : (3,) updated position
        v_n_out : (3,) updated velocity
        q_out   : (4,) updated quaternion
        """
        r_n = np.asarray(r_n, dtype=float).reshape(3)
        v_n = np.asarray(v_n, dtype=float).reshape(3)
        q = np.asarray(q, dtype=float).reshape(4)
        f_b = np.asarray(f_b, dtype=float).reshape(3)
        omega_ib_b = np.asarray(omega_ib_b, dtype=float).reshape(3)

        # -------------------------
        # 1) Attitude (quaternion)
        # -------------------------
        omega_e = 7.2921158e-5  # Earth rotation [rad/s]
        M, N = calM_N(r_n[0])

        omega_ie_n = np.array(
            [
                omega_e * np.cos(r_n[0]),
                0.0,
                -omega_e * np.sin(r_n[0]),
            ],
            dtype=float,
        )

        omega_en_n = np.array(
            [
                v_n[1] / (N + r_n[2]),
                -v_n[0] / (M + r_n[2]),
                -v_n[1] * np.tan(r_n[0]) / (N + r_n[2]),
            ],
            dtype=float,
        )

        dteta_ib_b = omega_ib_b * dt

        # Body-to-nav DCM from quaternion
        C_nb = CoordTransform.q2mat(q).T  # transpose to match MATLAB convention

        # Effective rotation of body wrt navigation frame over dt
        dteta_nb_b = dteta_ib_b - C_nb @ (omega_ie_n + omega_en_n) * dt
        dteta = np.sqrt(np.sum(dteta_nb_b**2))

        # Small-angle handling to avoid division by zero
        if dteta == 0.0:
            s = 0.0
            c = 0.0
        else:
            s = 2.0 * np.sin(dteta / 2.0) / dteta
            c = 2.0 * (np.cos(dteta / 2.0) - 1.0)

        # Quaternion update matrix
        temp_q = np.array(
            [
                [c, s * dteta_nb_b[2], -s * dteta_nb_b[1], s * dteta_nb_b[0]],
                [-s * dteta_nb_b[2], c, s * dteta_nb_b[0], s * dteta_nb_b[1]],
                [s * dteta_nb_b[1], -s * dteta_nb_b[0], c, s * dteta_nb_b[2]],
                [-s * dteta_nb_b[0], -s * dteta_nb_b[1], -s * dteta_nb_b[2], c],
            ],
            dtype=float,
        )

        q = q + 0.5 * (temp_q @ q)
        q = q / np.linalg.norm(q)

        # -------------------------
        # 2) Velocity in NED frame
        # -------------------------
        temp_v = np.array(
            [
                [1.0, dteta_nb_b[2] / 2.0, -dteta_nb_b[1] / 2.0],
                [-dteta_nb_b[2] / 2.0, 1.0, dteta_nb_b[0] / 2.0],
                [dteta_nb_b[1] / 2.0, -dteta_nb_b[0] / 2.0, 1.0],
            ],
            dtype=float,
        )

        dv_fb = f_b * dt
        dv_fn = C_nb.T @ (temp_v @ dv_fb)

        # Normal gravity model (height + latitude dependent)
        a1 = 9.7803267715
        a2 = 0.0052790414
        a3 = 0.0000232718
        a4 = -0.0000030876910891
        a5 = 0.0000000043977311
        a6 = 0.0000000000007211

        sin_lat = np.sin(r_n[0])
        gamma = (
            a1 * (1 + a2 * sin_lat**2 + a3 * sin_lat**4)
            + (a4 + a5 * sin_lat**2) * r_n[2]
            + a6 * r_n[2] ** 2
        )

        gamma_n = np.array([0.0, 0.0, gamma], dtype=float)

        dv_n = dv_fn - SkewSymmetry(2.0 * omega_ie_n + omega_en_n) @ v_n * dt + gamma_n * dt

        v_n_old = v_n.copy()
        v_n = v_n + dv_n

        # -------------------------
        # 3) Position update
        # -------------------------
        temp_r = np.array(
            [
                [1.0 / (M + r_n[2]), 0.0, 0.0],
                [0.0, 1.0 / ((N + r_n[2]) * np.cos(r_n[0])), 0.0],
                [0.0, 0.0, -1.0],
            ],
            dtype=float,
        )

        r_n = r_n + 0.5 * (temp_r @ (v_n_old + v_n) * dt)

        return r_n, v_n, q


