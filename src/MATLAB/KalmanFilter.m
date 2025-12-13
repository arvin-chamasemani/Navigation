classdef KalmanFilter
    %KALMANFILTER   Simple linear/extended Kalman filter class
    
    properties
        F   % State transition matrix (or Jacobian)
        Q   % Process noise covariance
        H   % Measurement matrix (or Jacobian)
        R   % Measurement noise covariance
        B   % Control input matrix (optional)
    end

    methods
        function obj = KalmanFilter(F, Q, H, R, B)
            % Constructor
            obj.F = F;
            obj.Q = Q;
            obj.H = H;
            obj.R = R;
            
            if nargin < 5
                obj.B = [];
            else
                obj.B = B;
            end
        end

        function [x_hat_new, P_new] = step(obj, x_hat, P, z, u)
            % STEP   Perform one predict + update cycle
            
            % ---------- Prediction ----------
            if isempty(obj.B)
                x_pred = obj.F * x_hat;
            else
                x_pred = obj.F * x_hat + obj.B * u;
            end
            
            P_pred = obj.F * P * obj.F' + obj.Q;
            
            % ---------- Update ----------
            y = z - obj.H * x_pred;                % innovation
            S = obj.H * P_pred * obj.H' + obj.R;   % innovation covariance
            K = P_pred * obj.H' / S;               % Kalman gain
            
            x_hat_new = x_pred + K * y;
            P_new     = (eye(size(P)) - K * obj.H) * P_pred;
        end
    end
    methods (Static)
        % function [x_hist, P_hist] = runBatchKF( ...
        %         t, r_cal, v_cal, q_cal, ...   % mechanization outputs
        %         r_meas,  ...          % "measurements" (e.g. r_n, v_n)
        %         f_b_IMU, ...     % IMU data
        %         Q_sens,var_GPS, ... % noise covariances
        %         P)                  % initial covariances 
        %     % RUNBATCHKF  Run error-state KF over the whole dataset
        %     %
        %     % Outputs:
        %     %   x_hist : [9 x N]  error-state history
        %     %   P_hist : [9 x 9 x N] covariance history
        %     %   r_corr : [3 x N]  corrected position (r_cal - δr)
        % 
        %     nx = 9;                     % [dr; dv; dtheta]
        %     N  = numel(t);
        % 
        %     % --- initial state and covariance ---
        %     x_hat = zeros(nx,1);
        % 
        %     % --- create KF object (F,Q,H will be overwritten every step) ---
        %     R     = diag(var_GPS);
        %     H = [eye(3), zeros(3),zeros(3)];    
        %     KF = KalmanFilter(eye(nx), eye(nx), H, R);
        % 
        %     % --- preallocate histories ---
        %     x_hist = zeros(nx, N);
        %     P_hist = zeros(nx, nx, N);
        % 
        % 
        %     % ================== MAIN KF LOOP ==================
        %     for k = 2:N
        %         dt = t(5) - t(4);
        % 
        %         % nav state from mechanization
        %         r_k = r_cal(:,k);
        %         v_k = v_cal(:,k);
        %         q_k = q_cal(:,k);
        % 
        %         % IMU data
        %         f_b_k      = f_b_IMU(:,k);
        % 
        %         % ---- build system matrices for this step ----
        %         F_k = KalmanFilter.F_cal(r_k, v_k, q_k, f_b_k, dt);
        %         Q_k = KalmanFilter.Q_cal(F_k, Q_sens, q_k, dt);
        % 
        %         KF.F = F_k;
        %         KF.Q = Q_k;
        %         KF.H = H;
        %         KF.R = R;
        % 
        %         % ---- measurement vector (position + velocity) ----
        %         z_k = r_meas(:,k)-r_cal(:,k);
        % 
        %         % ---- one KF step (returns ONLY x and P) ----
        %         [x_hat, P] = KF.step(x_hat, P, z_k, []);
        % 
        %         % store histories
        %         x_hist(:,k)     = x_hat;
        %         P_hist(:,:,k)   = P;
        % 
        %     end
    
        % end
    
        function F=F_cal(r_n,v_n,q,f_b,dt)
            omega_e=7.2921158e-5;
            [M,N]=Mechanization.calM_N(r_n(1));
            omega_ie_n = [omega_e*cos(r_n(1)), 0, -omega_e*sin(r_n(1))]';
            omega_en_n=[v_n(2)/(N+r_n(3)), -v_n(1)/(M+r_n(3)), -v_n(2)*tan(r_n(1))/(N+r_n(3))]';
            omega_in_n=omega_ie_n + omega_en_n;
            f_n=CoordTransform.q2mat(q)*f_b;

            F_rr= [0, 0, -v_n(1)/(M+r_n(3))^2;
                  v_n(2)*sin(r_n(1))/((N+r_n(3))*cos(r_n(1))^2),0,-v_n(2)/((N+r_n(3))^2*cos(r_n(1)));
                  0,0,0];
            F_rv= [1/(M+r_n(3)), 0, 0;
                   0, 1/(N+r_n(3))/cos(r_n(1)), 0;
                   0, 0, -1];
           
            F_vr= zeros(3,3);
            F_vr(1,1)= -2*v_n(2)*omega_e*cos(r_n(1))...
                - v_n(2)^2/((N+r_n(3))*cos(r_n(1))^2);
            F_vr(1,3)= -v_n(1)*v_n(3)/(M+r_n(3))^2;
            F_vr(2,1)= 2*omega_e*(v_n(1)*cos(r_n(1)) - v_n(3)*sin(r_n(1)))...
                + v_n(2)*v_n(1)/((N+r_n(3))*cos(r_n(1))^2);
            F_vr(2,3)= -v_n(2)*v_n(3)/(N+r_n(3))^2 - v_n(1)*v_n(2)*tan(r_n(1))/(N+r_n(3))^2;
            F_vr(3,1)= 2*v_n(2)*omega_e*sin(r_n(1));
            F_vr(3,3)= v_n(2)^2/(N+r_n(3))^2 - v_n(1)^2/(M+r_n(3))^2 ...
                - 2*KalmanFilter.gamma_cal(r_n(1),r_n(3))/(sqrt(M*N)+r_n(3));
            F_vv=zeros(3,3);
            F_vv(1,1) = v_n(3)/(M+r_n(3));
            F_vv(1,2) = -2*omega_e*sin(r_n(1)) - 2*v_n(2)*tan(r_n(1))/(N+r_n(3));
            F_vv(1,3) = v_n(1)/(M+r_n(3));
            F_vv(2,1) = 2*omega_e*sin(r_n(1)) + v_n(2)*tan(r_n(1))/(N+r_n(3)); 
            F_vv(2,2) = (v_n(3) + v_n(1)*tan(r_n(1)))/(N+r_n(3));
            F_vv(2,3) = 2*omega_e*cos(r_n(1))+v_n(2)/(N+r_n(3));
            F_vv(3,1) = -2*v_n(1)/(M+r_n(3));
            F_vv(3,2) = -2*omega_e*cos(r_n(1)) - 2*v_n(2)/(N+r_n(3));

            F_er =[-omega_e*sin(r_n(1)), 0, -v_n(2)/(N+r_n(3))^2;
                    0, 0, v_n(1)/(M+r_n(3))^2;
                    -omega_e*cos(r_n(1))-v_n(2)/(N+r_n(3))/cos(r_n(1))^2, 0, v_n(2)*tan(r_n(1))/(N+r_n(3))^2];
            F_ev= [0, 1/(N+r_n(3)),0;
                    -1/(M+r_n(3)), 0, 0;
                    0, -tan(r_n(1))/(N+r_n(3)), 0];

            F_0=[F_rr, F_rv, zeros(3,3);
                F_vr, F_vv, Mechanization.SkewSymmetry(f_n);
                F_er, F_ev, -1*Mechanization.SkewSymmetry(omega_in_n)];
        
            F=expm(F_0*dt);
            % F= eye(9)+F*dt;

        end
        function Q_k= Q_cal(F,Q,q,dt)
            C_bn=CoordTransform.q2mat(q);
            G=[zeros(3,6);
                C_bn, zeros(3,3);
                zeros(3,3), -C_bn];
            Q_k=F*G*Q*G'*F'*dt;
        end

        function gamma=gamma_cal(phi,h)
            a1 = 9.7803267715;a2 = 0.0052790414;a3 = 0.0000232718;
            a4 =-0.0000030876910891;a5 = 0.0000000043977311;a6 = 0.0000000000007211;

            gamma=a1*(1+a2*sin(phi)^2+a3*sin(phi)^4)+ (a4+a5*sin(phi)^2)*h+a6*h^2;
        end
    end
end