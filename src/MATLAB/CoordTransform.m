classdef CoordTransform
    % CoordTransform
    %   All routines assume ZYX Euler sequence (yaw->pitch->roll) and that
    %   DCMs map body-frame vectors to navigation-frame vectors:
    %       v_nav = C_nb * v_body
    %
    %   Quaternion format is q = [e1; e2; e3; eta] = [x; y; z; w]
    %   where eta = scalar part.
    methods(Static)
        %% --- Quaternion operations ---
        function q = qmul(q1, q2)
            %QMUL Multiply (compose) two quaternions.
            %
            %   q = CoordTransform.qmul(q1, q2)
            %
            %   Quaternion format: q = [x; y; z; w], w = scalar part.
            %
            %   This implements the Hamilton product:
            %       q = q1 ⊗ q2
            %
            %   If q1 and q2 represent rotations, the resulting rotation
            %   corresponds to "apply q2, then apply q1".
            %
            %   All inputs are normalized internally.

            % Ensure unit quaternions
            q1 = CoordTransform.normalizeQ(q1);
            q2 = CoordTransform.normalizeQ(q2);

            % Split into vector (v) and scalar (s) parts
            v1 = q1(1:3);  s1 = q1(4);
            v2 = q2(1:3);  s2 = q2(4);

            % Hamilton product
            v = s1*v2 + s2*v1 + cross(v1, v2);
            s = s1*s2 - dot(v1, v2);

            q = CoordTransform.normalizeQ([v; s]);
        end

        %% --- Utilities ---
        function q = normalizeQ(q)
            % Normalize a quaternion to unit length, with epsilon safety.
            n = norm(q);
            if n < eps
                q = [0;0;0;1];     % fall back to identity rotation
            else
                q = q / n;
            end
        end

        function x = clamp(x, lo, hi)
            % Clamp helper (scalar or array).
            x = min(max(x, lo), hi);
        end

        %% --- Euler (ZYX) <-> DCM ---
        % Convert Euler angles (ZYX: yaw-pitch-roll) to rotation matrix C_nb
        function C_nb = euler2mat(psi, teta, phi)
            % INPUTS:
            %   psi   - yaw   [rad] (about Z)
            %   teta  - pitch [rad] (about Y)
            %   phi   - roll  [rad] (about X)
            % OUTPUT:
            %   C_nb  - DCM mapping body -> nav: v_nav = C_nb * v_body

            cpsi  = cos(psi);   spsi  = sin(psi);
            cteta = cos(teta);  steta = sin(teta);
            cphi  = cos(phi);   sphi  = sin(phi);

            C_nb = [ ...
                cpsi*cteta,                      cpsi*steta*sphi - spsi*cphi,  cpsi*steta*cphi + spsi*sphi; ...
                spsi*cteta,                      spsi*steta*sphi + cpsi*cphi,  spsi*steta*cphi - cpsi*sphi; ...
                -steta,                          cteta*sphi,                   cteta*cphi                    ...
            ];
        end

        % Convert DCM (body->nav) to Euler angles (ZYX)
        function [psi, teta, phi] = mat2euler(C_nb)
            % INPUT:
            %   C_nb  - DCM mapping body -> nav
            % OUTPUT:
            %   psi   - yaw   [rad]
            %   teta  - pitch [rad]
            %   phi   - roll  [rad]

            % Numerical guard for asin
            t = CoordTransform.clamp(-C_nb(3,1), -1, 1); % -sin(theta)

            teta = asin(t);                     % pitch
            psi  = atan2(C_nb(2,1), C_nb(1,1)); % yaw
            phi  = atan2(C_nb(3,2), C_nb(3,3)); % roll
        end

        %% --- DCM <-> Quaternion (x,y,z,w with w=eta) ---
        function q = mat2q(C_nb)
            % Convert DCM (body->nav) to quaternion [x;y;z;w]
            % Reference: standard Tait-Bryan/ZYX mapping for passive DCM.
            tr = trace(C_nb);
            if tr > 0
                S  = sqrt(tr + 1.0) * 2;   % S = 4*w
                w  = 0.25 * S;
                x  = (C_nb(3,2) - C_nb(2,3)) / S;
                y  = (C_nb(1,3) - C_nb(3,1)) / S;
                z  = (C_nb(2,1) - C_nb(1,2)) / S;
            else
                % Find major diagonal term
                if (C_nb(1,1) > C_nb(2,2)) && (C_nb(1,1) > C_nb(3,3))
                    S  = sqrt(1.0 + C_nb(1,1) - C_nb(2,2) - C_nb(3,3)) * 2;
                    w  = (C_nb(3,2) - C_nb(2,3)) / S;
                    x  = 0.25 * S;
                    y  = (C_nb(1,2) + C_nb(2,1)) / S;
                    z  = (C_nb(1,3) + C_nb(3,1)) / S;
                elseif C_nb(2,2) > C_nb(3,3)
                    S  = sqrt(1.0 + C_nb(2,2) - C_nb(1,1) - C_nb(3,3)) * 2;
                    w  = (C_nb(1,3) - C_nb(3,1)) / S;
                    x  = (C_nb(1,2) + C_nb(2,1)) / S;
                    y  = 0.25 * S;
                    z  = (C_nb(2,3) + C_nb(3,2)) / S;
                else
                    S  = sqrt(1.0 + C_nb(3,3) - C_nb(1,1) - C_nb(2,2)) * 2;
                    w  = (C_nb(2,1) - C_nb(1,2)) / S;
                    x  = (C_nb(1,3) + C_nb(3,1)) / S;
                    y  = (C_nb(2,3) + C_nb(3,2)) / S;
                    z  = 0.25 * S;
                end
            end
            q = CoordTransform.normalizeQ([x; y; z; w]);
        end

        function C_nb = q2mat(q)
            % Convert quaternion [x;y;z;w] to DCM (body->nav)
            q = CoordTransform.normalizeQ(q);
            x = q(1); y = q(2); z = q(3); w = q(4);

            xx = x*x; yy = y*y; zz = z*z;
            xy = x*y; xz = x*z; yz = y*z;
            wx = w*x; wy = w*y; wz = w*z;

            C_nb = [ ...
                1 - 2*(yy + zz),   2*(xy - wz),       2*(xz + wy); ...
                2*(xy + wz),       1 - 2*(xx + zz),   2*(yz - wx); ...
                2*(xz - wy),       2*(yz + wx),       1 - 2*(xx + yy) ...
            ];
        end

        %% --- Quaternion <-> Euler (ZYX) ---
        % Preferred: returns yaw-pitch-roll (ZYX) in that order.
        function [psi, teta, phi] = q2eulerZYX(q)
            % OUTPUTS: psi (yaw), teta (pitch), phi (roll)
            q = CoordTransform.normalizeQ(q);
            x = q(1); y = q(2); z = q(3); w = q(4);

            % roll (phi)
            sinr_cosp = 2*(w*x + y*z);
            cosr_cosp = 1 - 2*(x*x + y*y);
            phi = atan2(sinr_cosp, cosr_cosp);

            % pitch (teta) with clamp
            sinp = 2*(w*y - z*x);
            sinp = CoordTransform.clamp(sinp, -1, 1);
            teta = asin(sinp);

            % yaw (psi)
            siny_cosp = 2*(w*z + x*y);
            cosy_cosp = 1 - 2*(y*y + z*z);
            psi = atan2(siny_cosp, cosy_cosp);
        end

        % Convenience: returns roll-pitch-yaw (RPY)
        function [phi, teta, psi] = q2rpy(q)
            [psi_, teta_, phi_] = CoordTransform.q2eulerZYX(q);
            phi  = phi_;
            teta = teta_;
            psi  = psi_;
        end

        % Back-compat signature (ZYX input order, returns [psi,teta,phi])
        function q = euler2q(psi, teta, phi)
            % INPUTS: psi (yaw), teta (pitch), phi (roll)
            cy = cos(psi/2);   sy = sin(psi/2);
            cp = cos(teta/2);  sp = sin(teta/2);
            cr = cos(phi/2);   sr = sin(phi/2);

            x = cy*cp*sr - sy*sp*cr;
            y = cy*sp*cr + sy*cp*sr;
            z = sy*cp*cr - cy*sp*sr;
            w = cy*cp*cr + sy*sp*sr;

            q = CoordTransform.normalizeQ([x; y; z; w]);
        end

        % Back-compat: returns [psi,teta,phi] (yaw,pitch,roll)
        function [psi, teta, phi] = q2euler(q)
            [psi, teta, phi] = CoordTransform.q2eulerZYX(q);
        end

        %% --- Demos & self-checks ---
        function demo()
            disp('--- Demo: Euler -> Matrix -> Quaternion -> Euler (ZYX) ---');
            psi = deg2rad(30); teta = deg2rad(60); phi = deg2rad(0);

            C = CoordTransform.euler2mat(psi, teta, phi);
            q = CoordTransform.mat2q(C);
            [psi2, teta2, phi2] = CoordTransform.q2eulerZYX(q);

            fprintf('Original  (deg): [yaw pitch roll] = [%g %g %g]\n', rad2deg([psi teta phi]));
            fprintf('Recovered (deg): [yaw pitch roll] = [%g %g %g]\n', rad2deg([psi2 teta2 phi2]));

            % quick consistency check
            C2 = CoordTransform.q2mat(q);
            err = norm(C - C2, 'fro');
            fprintf('||C - C(q)||_F = %.3e\n', err);
        end

        function demo2()
            disp('--- Demo2: Euler -> Quaternion -> Matrix -> Euler (ZYX) ---');
            psi = deg2rad(45); teta = deg2rad(15); phi = deg2rad(-30);

            q = CoordTransform.euler2q(psi, teta, phi);
            C = CoordTransform.q2mat(q);
            [psi2, teta2, phi2] = CoordTransform.mat2euler(C);

            fprintf('Original  (deg): [yaw pitch roll] = [%g %g %g]\n', rad2deg([psi teta phi]));
            fprintf('Recovered (deg): [yaw pitch roll] = [%g %g %g]\n', rad2deg([psi2 teta2 phi2]));
            err = norm(C - CoordTransform.euler2mat(psi2, teta2, phi2), 'fro');
            fprintf('||C - C(euler)||_F = %.3e\n', err);
        end
    end
end
