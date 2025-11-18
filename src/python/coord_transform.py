import numpy as np


class CoordTransform:
    """
    CoordTransform

    All routines assume ZYX Euler sequence (yaw -> pitch -> roll) and that
    DCMs map body-frame vectors to navigation-frame vectors:

        v_nav = C_nb @ v_body

    Quaternion format is q = [x, y, z, w] where w is the scalar part.
    """

    # --- Utilities ---

    @staticmethod
    def normalize_q(q):
        """
        Normalize a quaternion to unit length, with epsilon safety.
        """
        q = np.asarray(q, dtype=float).ravel()
        n = np.linalg.norm(q)
        if n < np.finfo(float).eps:
            return np.array([0.0, 0.0, 0.0, 1.0], dtype=float)
        return q / n

    @staticmethod
    def clamp(x, lo, hi):
        """
        Clamp helper (scalar or array).
        """
        x = np.asarray(x)
        return np.minimum(np.maximum(x, lo), hi)

    # --- Euler (ZYX) <-> DCM ---

    @staticmethod
    def euler2mat(psi, teta, phi):
        """
        Convert Euler angles (ZYX: yaw-pitch-roll) to rotation matrix C_nb.

        Inputs
        ------
        psi  : yaw   [rad] (about Z)
        teta : pitch [rad] (about Y)
        phi  : roll  [rad] (about X)

        Returns
        -------
        C_nb : (3,3) DCM mapping body -> nav: v_nav = C_nb @ v_body
        """
        cpsi = np.cos(psi)
        spsi = np.sin(psi)
        cteta = np.cos(teta)
        steta = np.sin(teta)
        cphi = np.cos(phi)
        sphi = np.sin(phi)

        C_nb = np.array([
            [cpsi * cteta,                    cpsi * steta * sphi - spsi * cphi,  cpsi * steta * cphi + spsi * sphi],
            [spsi * cteta,                    spsi * steta * sphi + cpsi * cphi,  spsi * steta * cphi - cpsi * sphi],
            [-steta,                          cteta * sphi,                       cteta * cphi]
        ], dtype=float)

        return C_nb

    @staticmethod
    def mat2euler(C_nb):
        """
        Convert DCM (body->nav) to Euler angles (ZYX).

        Input
        -----
        C_nb : (3,3) DCM mapping body -> nav

        Returns
        -------
        psi  : yaw   [rad]
        teta : pitch [rad]
        phi  : roll  [rad]
        """
        C_nb = np.asarray(C_nb, dtype=float).reshape(3, 3)

        # Numerical guard for asin
        t = CoordTransform.clamp(-C_nb[2, 0], -1.0, 1.0)  # -sin(theta)

        teta = np.arcsin(t)                     # pitch
        psi = np.arctan2(C_nb[1, 0], C_nb[0, 0])  # yaw
        phi = np.arctan2(C_nb[2, 1], C_nb[2, 2])  # roll

        return psi, teta, phi

    # --- DCM <-> Quaternion (x,y,z,w with w=eta) ---

    @staticmethod
    def mat2q(C_nb):
        """
        Convert DCM (body->nav) to quaternion [x, y, z, w].

        Reference: standard Tait-Bryan/ZYX mapping for passive DCM.
        """
        C_nb = np.asarray(C_nb, dtype=float).reshape(3, 3)
        tr = np.trace(C_nb)

        if tr > 0.0:
            S = np.sqrt(tr + 1.0) * 2.0  # S = 4*w
            w = 0.25 * S
            x = (C_nb[2, 1] - C_nb[1, 2]) / S
            y = (C_nb[0, 2] - C_nb[2, 0]) / S
            z = (C_nb[1, 0] - C_nb[0, 1]) / S
        else:
            if (C_nb[0, 0] > C_nb[1, 1]) and (C_nb[0, 0] > C_nb[2, 2]):
                S = np.sqrt(1.0 + C_nb[0, 0] - C_nb[1, 1] - C_nb[2, 2]) * 2.0
                w = (C_nb[2, 1] - C_nb[1, 2]) / S
                x = 0.25 * S
                y = (C_nb[0, 1] + C_nb[1, 0]) / S
                z = (C_nb[0, 2] + C_nb[2, 0]) / S
            elif C_nb[1, 1] > C_nb[2, 2]:
                S = np.sqrt(1.0 + C_nb[1, 1] - C_nb[0, 0] - C_nb[2, 2]) * 2.0
                w = (C_nb[0, 2] - C_nb[2, 0]) / S
                x = (C_nb[0, 1] + C_nb[1, 0]) / S
                y = 0.25 * S
                z = (C_nb[1, 2] + C_nb[2, 1]) / S
            else:
                S = np.sqrt(1.0 + C_nb[2, 2] - C_nb[0, 0] - C_nb[1, 1]) * 2.0
                w = (C_nb[1, 0] - C_nb[0, 1]) / S
                x = (C_nb[0, 2] + C_nb[2, 0]) / S
                y = (C_nb[1, 2] + C_nb[2, 1]) / S
                z = 0.25 * S

        q = np.array([x, y, z, w], dtype=float)
        return CoordTransform.normalize_q(q)

    @staticmethod
    def q2mat(q):
        """
        Convert quaternion [x, y, z, w] to DCM (body->nav).
        """
        q = CoordTransform.normalize_q(q)
        x, y, z, w = q

        xx = x * x
        yy = y * y
        zz = z * z
        xy = x * y
        xz = x * z
        yz = y * z
        wx = w * x
        wy = w * y
        wz = w * z

        C_nb = np.array([
            [1.0 - 2.0 * (yy + zz),   2.0 * (xy - wz),         2.0 * (xz + wy)],
            [2.0 * (xy + wz),         1.0 - 2.0 * (xx + zz),   2.0 * (yz - wx)],
            [2.0 * (xz - wy),         2.0 * (yz + wx),         1.0 - 2.0 * (xx + yy)]
        ], dtype=float)

        return C_nb

    # --- Quaternion <-> Euler (ZYX) ---

    @staticmethod
    def q2eulerZYX(q):
        """
        Preferred: returns yaw-pitch-roll (ZYX) in that order.

        Outputs
        -------
        psi  : yaw   [rad]
        teta : pitch [rad]
        phi  : roll  [rad]
        """
        q = CoordTransform.normalize_q(q)
        x, y, z, w = q

        # roll (phi)
        sinr_cosp = 2.0 * (w * x + y * z)
        cosr_cosp = 1.0 - 2.0 * (x * x + y * y)
        phi = np.arctan2(sinr_cosp, cosr_cosp)

        # pitch (teta) with clamp
        sinp = 2.0 * (w * y - z * x)
        sinp = CoordTransform.clamp(sinp, -1.0, 1.0)
        teta = np.arcsin(sinp)

        # yaw (psi)
        siny_cosp = 2.0 * (w * z + x * y)
        cosy_cosp = 1.0 - 2.0 * (y * y + z * z)
        psi = np.arctan2(siny_cosp, cosy_cosp)

        return psi, teta, phi

    @staticmethod
    def q2rpy(q):
        """
        Convenience: returns roll-pitch-yaw (RPY).
        """
        psi_, teta_, phi_ = CoordTransform.q2eulerZYX(q)
        phi = phi_
        teta = teta_
        psi = psi_
        return phi, teta, psi

    @staticmethod
    def euler2q(psi, teta, phi):
        """
        Back-compat signature (ZYX input order, returns [x, y, z, w]).

        Inputs
        ------
        psi  : yaw   [rad]
        teta : pitch [rad]
        phi  : roll  [rad]
        """
        cy = np.cos(psi / 2.0)
        sy = np.sin(psi / 2.0)
        cp = np.cos(teta / 2.0)
        sp = np.sin(teta / 2.0)
        cr = np.cos(phi / 2.0)
        sr = np.sin(phi / 2.0)

        x = cy * cp * sr - sy * sp * cr
        y = cy * sp * cr + sy * cp * sr
        z = sy * cp * cr - cy * sp * sr
        w = cy * cp * cr + sy * sp * sr

        q = np.array([x, y, z, w], dtype=float)
        return CoordTransform.normalize_q(q)

    @staticmethod
    def q2euler(q):
        """
        Back-compat: returns [psi, teta, phi] (yaw, pitch, roll).
        """
        return CoordTransform.q2eulerZYX(q)

    # --- Demos & self-checks ---

    @staticmethod
    def demo():
        """
        Demo: Euler -> Matrix -> Quaternion -> Euler (ZYX).
        """
        print('--- Demo: Euler -> Matrix -> Quaternion -> Euler (ZYX) ---')
        psi = np.deg2rad(30.0)
        teta = np.deg2rad(60.0)
        phi = np.deg2rad(0.0)

        C = CoordTransform.euler2mat(psi, teta, phi)
        q = CoordTransform.mat2q(C)
        psi2, teta2, phi2 = CoordTransform.q2eulerZYX(q)

        original = np.rad2deg([psi, teta, phi])
        recovered = np.rad2deg([psi2, teta2, phi2])

        print(f'Original  (deg): [yaw pitch roll] = {original}')
        print(f'Recovered (deg): [yaw pitch roll] = {recovered}')

        C2 = CoordTransform.q2mat(q)
        err = np.linalg.norm(C - C2, ord='fro')
        print(f'||C - C(q)||_F = {err:.3e}')

    @staticmethod
    def demo2():
        """
        Demo2: Euler -> Quaternion -> Matrix -> Euler (ZYX).
        """
        print('--- Demo2: Euler -> Quaternion -> Matrix -> Euler (ZYX) ---')
        psi = np.deg2rad(45.0)
        teta = np.deg2rad(15.0)
        phi = np.deg2rad(-30.0)

        q = CoordTransform.euler2q(psi, teta, phi)
        C = CoordTransform.q2mat(q)
        psi2, teta2, phi2 = CoordTransform.mat2euler(C)

        original = np.rad2deg([psi, teta, phi])
        recovered = np.rad2deg([psi2, teta2, phi2])

        print(f'Original  (deg): [yaw pitch roll] = {original}')
        print(f'Recovered (deg): [yaw pitch roll] = {recovered}')

        C2 = CoordTransform.euler2mat(psi2, teta2, phi2)
        err = np.linalg.norm(C - C2, ord='fro')
        print(f'||C - C(euler)||_F = {err:.3e}')


if __name__ == "__main__":
    # Quick self-test if you run this file directly
    CoordTransform.demo()
    CoordTransform.demo2()
