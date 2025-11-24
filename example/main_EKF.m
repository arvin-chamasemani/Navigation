clc
clear
close all
% programmer: Arvin Chamasemani

% ensure script runs from its own folder
cd(fileparts(mfilename('fullpath')))

% loading data
load("../data/path2sensor_output.mat");

% main class
[r_cal, v_cal, q_cal, roll_cal, pitch_cal, yaw_cal] = ...
    Mechanization.runBatch(t, r_n, v_n, roll, pitch, yaw, f_b_IMU, w_b_ib_IMU);

% Kalman filter calculations
P = diag([1,1,1,1e-3,1e-3,1e-3,1e-9,1e-9,1e-9]); % initial covariances

[x_hist, P_hist] = KalmanFilter.runBatchKF( ...
    t, r_cal, v_cal, q_cal, ...     % mechanization outputs
    r_n,  ...                   % "measurements"
    f_b_IMU, ...        % IMU data
    Q, var_GPS, ...      % noise covariances
    P);                 % initial covariances   
% updating states
r_ekf = zeros(3,numel(t));
v_ekf = zeros(3,numel(t));

for i= 1:numel(t)
    r_ekf(:,i)=r_cal(:,i)+x_hist(1:3,i);
    v_ekf(:,i)=v_cal(:,i)+x_hist(4:6,i);
    C_bn_corr= (eye(3)+Mechanization.SkewSymmetry(x_hist(7:9,i)))*CoordTransform.q2mat(q_cal(:,i));
    [yaw_ekf(i),pitch_ekf(i),roll_ekf(i)]=CoordTransform.mat2euler(C_bn_corr);
end

% ================== PLOTTING ==================
plotting_ekf( t, ...
          r_n, v_n, roll, pitch, yaw, ...     % reference
          r_cal, v_cal, roll_cal, pitch_cal, yaw_cal,... % INS-based output
          r_ekf, v_ekf, roll_ekf, pitch_ekf, yaw_ekf );   % EKF output
