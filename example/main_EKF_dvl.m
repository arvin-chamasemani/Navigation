%%
clc
clear
close all
% programmer: Arvin Chamasemani

cd(fileparts(mfilename('fullpath')))

% load data
load("../data/path2sensor_output.mat");

N = numel(t);

% --- initial mechanized navigation state ---
r_n_cal(:,1) = r_n(:,1);
v_n_cal(:,1) = v_n(:,1);
q_cal(:,1)   = CoordTransform.euler2q(yaw(1), pitch(1), roll(1));

% --- EKF initialization ---
nx = 9;     % δr, δv, δθ
x_hat = zeros(nx,1);

[M,N_n]=Mechanization.calM_N(r_n(1));

P = diag([1e-15,1e-15,1, 1e-3,1e-3,1e-3, 1e-9,1e-9,1e-9]);

% --- storage ---
x_hist   = zeros(nx, N);
P_hist   = zeros(nx, nx, N);

r_feedback_ekf   = zeros(3,N);
v_feedback_ekf   = zeros(3,N);
roll_feedback_ekf  = zeros(1,N);
pitch_feedback_ekf = zeros(1,N);
yaw_feedback_ekf   = zeros(1,N);

% initialize feedback logs
r_feedback_ekf(:,1) = r_n_cal(:,1);
v_feedback_ekf(:,1) = v_n_cal(:,1);
[yaw_feedback_ekf(1), pitch_feedback_ekf(1), roll_feedback_ekf(1)] = ...
    CoordTransform.q2euler(q_cal(:,1));

%% ================== FEEDBACK EKF LOOP ==================
for i = 2:N
    
    dt = t(5) - t(4);

    % --- 1) Mechanization ---
    [r_n_cal(:,i), v_n_cal(:,i), q_cal(:,i)] = Mechanization.dynamics( ...
        r_n_cal(:,i-1), v_n_cal(:,i-1), q_cal(:,i-1), ...
        f_b_IMU(:,i-1), w_b_ib_IMU(:,i-1), dt);
    C_bn = CoordTransform.q2mat(q_cal(:,i));

    r_k = r_n_cal(:,i);
    v_k = v_n_cal(:,i);
    q_k = q_cal(:,i);

    % --- IMU force for Jacobian ---
    f_b_k = f_b_IMU(:,i-1);

    % --- 2) Build step matrices ---
    % R = diag([var_GPS;var_DVL']); 
      R = diag([1e-3;1e-3;1e-3;1e-3;1e-3;1e-3]);                  % GPS output covariance
    % GPS output covariance
    H = [eye(3), zeros(3), zeros(3);
    zeros(3), eye(3), -Mechanization.SkewSymmetry(C_bn*v_DVL(:,i))];   % position-only update
    % tmp=[(M+r_n(3)), 0, 0 ;
    %     0, (N_n+r_n(3))*cos(r_n(1)),0;
    %     0, 0, 1];
    % H = [tmp , zeros(3,6);
    %     zeros(3), eye(3), -Mechanization.SkewSymmetry(C_bn*v_DVL(:,i))];   % position-only update

    F_k = KalmanFilter.F_cal(r_k, v_k, q_k, f_b_k, dt);
    Q_k = KalmanFilter.Q_cal(F_k, Q, q_k, dt);

    KF = KalmanFilter(eye(nx), eye(nx), H, R);

    KF.F = F_k;
    KF.Q = Q_k;
    KF.H = H;
    KF.R = R;
    
    % --- 3) Measurement residual (position only) ---
    z_k = [(r_n_cal(:,i)-GPS(:,i));
        (v_n_cal(:,i)-C_bn*v_DVL(:,i))];
    % z_k = [(M+r_n(3))*(r_n_cal(1,i)-GPS(1,i));
    %     (N_n+r_n(3))*cos(r_n(1))*(r_n_cal(2,i)-GPS(2,i));
    %     r_n_cal(3,i)-GPS(3,i);
    %     (v_n_cal(:,i)-C_bn*v_DVL(:,i))];

    % --- 4) EKF step ---
    [x_hat, P] = KF.step(x_hat, P, z_k, []);

    x_hist(:,i)   = x_hat;
    P_hist(:,:,i) = P;

    % ================== 5) FEEDBACK / STATE INJECTION ==================

    % position + velocity correction
    r_n_cal(:,i) = r_n_cal(:,i) - x_hat(1:3);
    v_n_cal(:,i) = v_n_cal(:,i) - x_hat(4:6);


    % attitude correction
    % C_corr = (eye(3) + Mechanization.SkewSymmetry(x_hat(7:9))) * C_bn;
    q_cal(:,i) = CoordTransform.qmul([0.5*x_hat(7:9);1],q_cal(:,i));
    [yaw_feedback_ekf(i), pitch_feedback_ekf(i), roll_feedback_ekf(i)] = ...
        CoordTransform.q2euler(q_cal(:,i));

    % q_cal(:,i) = CoordTransform.mat2q( C_corr);

    % store corrected feedback trajectory
    r_feedback_ekf(:,i) = r_n_cal(:,i);
    v_feedback_ekf(:,i) = v_n_cal(:,i);

    % IMPORTANT: Reset error state after feedback
    x_hat = zeros(nx,1);

end



%% === Run batch INS  ===
[r_cal, v_cal, q_cal, yaw_cal, pitch_cal, roll_cal] = ...
    Mechanization.runBatch(t, r_n, v_n, yaw, pitch, roll, f_b_IMU, w_b_ib_IMU);

%% === Run feedforward EKF  ===

% Kalman filter calculations
P = diag([1e-15,1e-15,1,1e-3,1e-3,1e-3,1e-9,1e-9,1e-9]); % initial covariances

% RUNBATCHKF  Run error-state KF over the whole dataset
%
% Outputs:
%   x_hist : [9 x N]  error-state history
%   P_hist : [9 x 9 x N] covariance history
%   r_corr : [3 x N]  corrected position (r_cal - δr)

nx = 9;                     % [dr; dv; dtheta]
N  = numel(t);

% --- initial state and covariance ---
x_hat = zeros(nx,1);

% --- create KF object (F,Q,H will be overwritten every step) ---
R     = diag(var_GPS);
H = [eye(3), zeros(3),zeros(3)];    
KF = KalmanFilter(eye(nx), eye(nx), H, R);

% --- preallocate histories ---
x_hist = zeros(nx, N);
P_hist = zeros(nx, nx, N);


% ================== MAIN KF LOOP ==================
for k = 2:N
    dt = t(5) - t(4);

    % nav state from mechanization
    r_k = r_cal(:,k);
    v_k = v_cal(:,k);
    q_k = q_cal(:,k);

    % IMU data
    f_b_k      = f_b_IMU(:,k-1);

    % ---- build system matrices for this step ----
    F_k = KalmanFilter.F_cal(r_k, v_k, q_k, f_b_k, dt);
    Q_k = KalmanFilter.Q_cal(F_k, Q, q_k, dt);

    KF.F = F_k;
    KF.Q = Q_k;
    KF.H = H;
    KF.R = R;

    % ---- measurement vector (position + velocity) ----
    z_k = r_cal(:,k)-GPS(:,k);

    % ---- one KF step (returns ONLY x and P) ----
    [x_hat, P] = KF.step(x_hat, P, z_k, []);

    % store histories
    x_hist(:,k)     = x_hat;
    P_hist(:,:,k)   = P;

end


r_ekf = zeros(3,numel(t));
v_ekf = zeros(3,numel(t));

for i= 1:numel(t)
    r_ekf(:,i) = r_cal(:,i)-x_hist(1:3,i);
    v_ekf(:,i) = v_cal(:,i)-x_hist(4:6,i);

   delta_theta = x_hist(7:9,i);                             % attitude error at time i
dq = CoordTransform.normalizeQ([0.5*delta_theta; 1]);    % small-angle quaternion
q_corr = CoordTransform.qmul(dq, q_cal(:,i));            % corrected attitude
[yaw_ekf(i), pitch_ekf(i), roll_ekf(i)] = ...
    CoordTransform.q2euler(q_corr);
end

%% === Plotting ===
plotting_ekf( t, ...
          r_n, v_n, roll, pitch, yaw, ...      
          r_cal, v_cal, roll_cal, pitch_cal, yaw_cal, ...
          r_ekf,v_ekf,roll_ekf,pitch_ekf,yaw_ekf,...  
          r_feedback_ekf, v_feedback_ekf,roll_feedback_ekf, pitch_feedback_ekf, yaw_feedback_ekf);
