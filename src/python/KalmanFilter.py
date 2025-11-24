import numpy as np
from scipy.linalg import expm

from mechanization import Mechanization_OOP      # must provide calM_N, skew_symmetric
from coord_transform import CoordTransform   # must provide q2mat(q)
from SkewSymmetry import SkewSymmetry
from calM_N import calM_N


class KalmanFilter:
    """
    Simple linear / error-state Kalman filter.

    Attributes
    ----------
    F : np.ndarray
        State transition matrix (or Jacobian).
    Q : np.ndarray
        Process noise covariance.
    H : np.ndarray
        Measurement matrix (or Jacobian).
    R : np.ndarray
        Measurement noise covariance.
    B : np.ndarray or None
        Control input matrix (optional).
    """

    def __init__(self, F, Q, H, R, B=None):
        self.F = F
        self.Q = Q
        self.H = H
        self.R = R
        self.B = B

    def step(self, x_hat, P, z, u=None):
        """
        Perform one predict + update cycle.

        Parameters
        ----------
        x_hat : np.ndarray, shape (nx,)
            Current state estimate.
        P : np.ndarray, shape (nx, nx)
            Current covariance.
        z : np.ndarray
            Measurement vector.
        u : np.ndarray or None
            Control input (optional).

        Returns
        -------
        x_hat_new : np.ndarray
            Updated state estimate.
        P_new : np.ndarray
            Updated covariance.
        """
        # ---------- Prediction ----------
        if self.B is None or u is None:
            x_pred = self.F @ x_hat
        else:
            x_pred = self.F @ x_hat + self.B @ u

        P_pred = self.F @ P @ self.F.T + self.Q

        # ---------- Update ----------
        y = z - (self.H @ x_pred)                 # innovation
        S = self.H @ P_pred @ self.H.T + self.R   # innovation covariance
        K = P_pred @ self.H.T @ np.linalg.inv(S)  # Kalman gain

        x_hat_new = x_pred + K @ y
        I = np.eye(P.shape[0])
        P_new = (I - K @ self.H) @ P_pred

        return x_hat_new, P_new

    # ==================================================================
    # Static helper routines, translated from your MATLAB code
    # ==================================================================

    @staticmethod
    def run_batch_kf(
        t,
        r_cal, v_cal, q_cal,      # mechanization outputs (3xN, 3xN, 4xN or 4xN)
        r_meas,                  # "measurements" (e.g., r_n) (3xN)
        f_b_IMU,                 # IMU specific force (3xN)
        Q_sens,                  # 6x6 sensor noise covariance
        var_GPS,                 # length-3 array for position measurement noise variances
        P                        # initial covariance (9x9)
    ):
        """
        Run error-state KF over the whole dataset (position-only measurement here).

        Parameters
        ----------
        t : np.ndarray, shape (N,)
            Time vector.
        r_cal, v_cal : np.ndarray, shape (3, N)
            Mechanization position and velocity.
        q_cal : np.ndarray, shape (4, N)
            Mechanization attitude (quaternions).
        r_meas : np.ndarray, shape (3, N)
            Measured position (e.g., from GPS).
        f_b_IMU : np.ndarray, shape (3, N)
            Body-frame specific force from IMU.
        Q_sens : np.ndarray, shape (6, 6)
            Sensor noise covariance.
        var_GPS : np.ndarray, shape (3,)
            Position measurement noise variances.
        P : np.ndarray, shape (9, 9)
            Initial state covariance.

        Returns
        -------
        x_hist : np.ndarray, shape (9, N)
            Error-state history [dr; dv; dtheta].
        P_hist : np.ndarray, shape (9, 9, N)
            Covariance history.
        """
        nx = 9                  # [dr; dv; dtheta]
        N = t.size

        # Initial state (zero error)
        x_hat = np.zeros(nx)

        # Initial measurement model: position only
        R = np.diag(var_GPS)
        H = np.hstack([
            np.eye(3),
            np.zeros((3, 3)),
            np.zeros((3, 3))
        ])   # [ I_3  0  0 ]

        # Create KF object (F and Q overwritten each step)
        KF = KalmanFilter(F=np.eye(nx), Q=np.eye(nx), H=H, R=R)

        # Preallocate histories
        x_hist = np.zeros((nx, N))
        P_hist = np.zeros((nx, nx, N))

        x_hist[:, 0] = x_hat
        P_hist[:, :, 0] = P

        # ================== MAIN KF LOOP ==================
        for k in range(1, N):
            dt = float(t[k] - t[k - 1])

            # nav state from mechanization
            r_k = r_cal[:, k]
            v_k = v_cal[:, k]
            q_k = q_cal[:, k]

            # IMU data
            f_b_k = f_b_IMU[:, k]

            # ---- build system matrices for this step ----
            F_k = KalmanFilter.F_cal(r_k, v_k, q_k, f_b_k, dt)
            Q_k = KalmanFilter.Q_cal(F_k, Q_sens, q_k, dt)

            KF.F = F_k
            KF.Q = Q_k
            KF.H = H
            KF.R = R

            # ---- measurement vector (position error) ----
            z_k = r_meas[:, k] - r_cal[:, k]

            # ---- one KF step ----
            x_hat, P = KF.step(x_hat, P, z_k, u=None)

            # store histories
            x_hist[:, k] = x_hat
            P_hist[:, :, k] = P

        return x_hist, P_hist

    @staticmethod
    def F_cal(r_n, v_n, q, f_b, dt):
        """
        Build discrete-time state transition matrix F_k (9x9) for one step.
        r_n, v_n : 3-vector (lat, lon, h) and (vN, vE, vD) in NED.
        q : 4-vector quaternion (b->n).
        f_b : 3-vector specific force in body frame.
        """
        omega_e = 7.2921158e-5
        phi = r_n[0]
        h = r_n[2]

        M, N_rad = calM_N(phi)  # must return (M, N)

        omega_ie_n = np.array([
            omega_e * np.cos(phi),
            0.0,
            -omega_e * np.sin(phi)
        ])

        omega_en_n = np.array([
            v_n[1] / (N_rad + h),
            -v_n[0] / (M + h),
            -v_n[1] * np.tan(phi) / (N_rad + h)
        ])

        omega_in_n = omega_ie_n + omega_en_n

        C_nb = CoordTransform.q2mat(q)  # 3x3
        f_n = C_nb @ f_b

        # --- Sub-blocks (3x3) ---
        # F_rr
        F_rr = np.zeros((3, 3))
        F_rr[0, 2] = -v_n[0] / (M + h) ** 2
        F_rr[1, 0] = v_n[1] * np.sin(phi) / ((N_rad + h) * np.cos(phi) ** 2)
        F_rr[1, 2] = -v_n[1] / ((N_rad + h) ** 2 * np.cos(phi))
        # F_rr[2,:] = 0

        # F_rv
        F_rv = np.zeros((3, 3))
        F_rv[0, 0] = 1.0 / (M + h)
        F_rv[1, 1] = 1.0 / ((N_rad + h) * np.cos(phi))
        F_rv[2, 2] = -1.0

        # F_vr
        F_vr = np.zeros((3, 3))
        F_vr[0, 0] = (
            -2.0 * v_n[1] * omega_e * np.cos(phi)
            - v_n[1] ** 2 / ((N_rad + h) * np.cos(phi) ** 2)
        )
        F_vr[0, 2] = -v_n[0] * v_n[2] / (M + h) ** 2

        F_vr[1, 0] = (
            2.0 * omega_e * (v_n[0] * np.cos(phi) - v_n[2] * np.sin(phi))
            + v_n[1] * v_n[0] / ((N_rad + h) * np.cos(phi) ** 2)
        )
        F_vr[1, 2] = (
            -v_n[1] * v_n[2] / (N_rad + h) ** 2
            - v_n[0] * v_n[1] * np.tan(phi) / (N_rad + h) ** 2
        )

        F_vr[2, 0] = 2.0 * v_n[1] * omega_e * np.sin(phi)
        F_vr[2, 2] = (
            v_n[1] ** 2 / (N_rad + h) ** 2
            - v_n[0] ** 2 / (M + h) ** 2
            - 2.0 * KalmanFilter.gamma_cal(phi, h) / (np.sqrt(M * N_rad) + h)
        )

        # F_vv
        F_vv = np.zeros((3, 3))
        F_vv[0, 0] = v_n[2] / (M + h)
        F_vv[0, 1] = -2.0 * omega_e * np.sin(phi) - 2.0 * v_n[1] * np.tan(phi) / (N_rad + h)
        F_vv[0, 2] = v_n[0] / (M + h)

        F_vv[1, 0] = 2.0 * omega_e * np.sin(phi) + v_n[1] * np.tan(phi) / (N_rad + h)
        F_vv[1, 1] = (v_n[2] + v_n[0] * np.tan(phi)) / (N_rad + h)
        F_vv[1, 2] = 2.0 * omega_e * np.cos(phi) + v_n[1] / (N_rad + h)

        F_vv[2, 0] = -2.0 * v_n[0] / (M + h)
        F_vv[2, 1] = -2.0 * omega_e * np.cos(phi) - 2.0 * v_n[1] / (N_rad + h)

        # F_er
        F_er = np.zeros((3, 3))
        F_er[0, 0] = -omega_e * np.sin(phi)
        F_er[0, 2] = -v_n[1] / (N_rad + h) ** 2

        F_er[1, 2] = v_n[0] / (M + h) ** 2

        F_er[2, 0] = -omega_e * np.cos(phi) - v_n[1] / ((N_rad + h) * np.cos(phi) ** 2)
        F_er[2, 2] = v_n[1] * np.tan(phi) / (N_rad + h) ** 2

        # F_ev
        F_ev = np.zeros((3, 3))
        F_ev[0, 1] = 1.0 / (N_rad + h)
        F_ev[1, 0] = -1.0 / (M + h)
        F_ev[2, 1] = -np.tan(phi) / (N_rad + h)

        # Assemble continuous-time F (9x9)
        F_cont = np.block([
            [F_rr,              F_rv,                       np.zeros((3, 3))],
            [F_vr,              F_vv,                       SkewSymmetry(f_n)],
            [F_er,              F_ev,                       SkewSymmetry(omega_in_n)],
        ])

        # Discretize
        F_disc = expm(F_cont * dt)
        # Or first-order: F_disc = np.eye(9) + F_cont * dt

        return F_disc

    @staticmethod
    def Q_cal(F, Q, q, dt):
        """
        Discrete process noise covariance Q_k.
        Q is the continuous sensor noise covariance (6x6).
        """
        C_nb = CoordTransform.q2mat(q)
        G = np.block([
            [np.zeros((3, 3)), np.zeros((3, 3))],
            [C_nb,             np.zeros((3, 3))],
            [np.zeros((3, 3)), -C_nb]
        ])   # shape (9, 6)

        Q_k = F @ G @ Q @ G.T @ F.T * dt
        return Q_k

    @staticmethod
    def gamma_cal(phi, h):
        """
        Gravity model as in your MATLAB code.
        """
        a1 = 9.7803267715
        a2 = 0.0052790414
        a3 = 0.0000232718
        a4 = -0.0000030876910891
        a5 = 0.0000000043977311
        a6 = 0.0000000000007211

        sin_phi = np.sin(phi)
        gamma = (
            a1 * (1.0 + a2 * sin_phi ** 2 + a3 * sin_phi ** 4)
            + (a4 + a5 * sin_phi ** 2) * h
            + a6 * h ** 2
        )
        return gamma
