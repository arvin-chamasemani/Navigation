function plotting_ekf(t, r_n, v_n, roll, pitch, yaw, ...
                  r_n_cal, v_n_cal, roll_cal, pitch_cal, yaw_cal, ...
                  r_n_ekf, v_n_ekf, roll_ekf, pitch_ekf, yaw_ekf, ...
                  r_fb_ekf, v_fb, roll_fb, pitch_fb, yaw_fb)
%PLOTTING  Compare reference vs INS vs EKF vs FB-EKF (velocity, position, attitude)
%          and plot errors (INS, EKF, FB-EKF) on same axes.
%          RMSEs are summarized in a separate figure as a table.

%% ===================== VELOCITY (N,E,D) ============================
% Errors: (estimate - reference)
vel_err_mech = v_n_cal - v_n;   % INS - ref
vel_err_ekf  = v_n_ekf - v_n;   % EKF - ref
vel_err_fb   = v_fb    - v_n;   % FB-EKF - ref

figure('Name','Velocity (NED) - Comparison & Error','Color','w');

% --- North velocity ---
subplot(3,2,1);
plot(t, v_n_cal(1,:), 'b',  'LineWidth', 1.2); hold on;
plot(t, v_n_ekf(1,:), 'g--','LineWidth', 1.2);
plot(t, v_fb(1,:),    'c-.','LineWidth', 1.2);
plot(t, v_n(1,:),     'r:', 'LineWidth', 1.2);
grid on;
legend('INS','EKF','FB EKF','Ref','Location','best');
title('North Velocity');
ylabel('v_N [m/s]');
xlabel('Time [s]');

subplot(3,2,2);
plot(t, vel_err_mech(1,:), 'k', 'LineWidth', 1.2); hold on;
plot(t, vel_err_ekf(1,:),  'm', 'LineWidth', 1.2);
plot(t, vel_err_fb(1,:),   'c--','LineWidth', 1.2);
grid on;
legend('INS Error','EKF Error','FB Error','Location','best');
title('North Velocity Error');
ylabel('Error [m/s]');
xlabel('Time [s]');
ylim(autoErrorYLim([vel_err_mech(1,:), vel_err_ekf(1,:), vel_err_fb(1,:)]));

% --- East velocity ---
subplot(3,2,3);
plot(t, v_n_cal(2,:), 'b',  'LineWidth', 1.2); hold on;
plot(t, v_n_ekf(2,:), 'g--','LineWidth', 1.2);
plot(t, v_fb(2,:),    'c-.','LineWidth', 1.2);
plot(t, v_n(2,:),     'r:', 'LineWidth', 1.2);
grid on;
legend('INS','EKF','FB EKF','Ref','Location','best');
title('East Velocity');
ylabel('v_E [m/s]');
xlabel('Time [s]');

subplot(3,2,4);
plot(t, vel_err_mech(2,:), 'k', 'LineWidth', 1.2); hold on;
plot(t, vel_err_ekf(2,:),  'm', 'LineWidth', 1.2);
plot(t, vel_err_fb(2,:),   'c--','LineWidth', 1.2);
grid on;
legend('INS Error','EKF Error','FB Error','Location','best');
title('East Velocity Error');
ylabel('Error [m/s]');
xlabel('Time [s]');
ylim(autoErrorYLim([vel_err_mech(2,:), vel_err_ekf(2,:), vel_err_fb(2,:)]));

% --- Down velocity ---
subplot(3,2,5);
plot(t, v_n_cal(3,:), 'b',  'LineWidth', 1.2); hold on;
plot(t, v_n_ekf(3,:), 'g--','LineWidth', 1.2);
plot(t, v_fb(3,:),    'c-.','LineWidth', 1.2);
plot(t, v_n(3,:),     'r:', 'LineWidth', 1.2);
grid on;
legend('INS','EKF','FB EKF','Ref','Location','best');
title('Down Velocity');
ylabel('v_D [m/s]');
xlabel('Time [s]');

subplot(3,2,6);
plot(t, vel_err_mech(3,:), 'k', 'LineWidth', 1.2); hold on;
plot(t, vel_err_ekf(3,:),  'm', 'LineWidth', 1.2);
plot(t, vel_err_fb(3,:),   'c--','LineWidth', 1.2);
grid on;
legend('INS Error','EKF Error','FB Error','Location','best');
title('Down Velocity Error');
ylabel('Error [m/s]');
xlabel('Time [s]');
ylim(autoErrorYLim([vel_err_mech(3,:), vel_err_ekf(3,:), vel_err_fb(3,:)]));

%% ===================== POSITION (Lat, Lon, H) ======================

% Convert lat/lon to degrees
lat_ref = rad2deg(r_n(1,:));
lon_ref = rad2deg(r_n(2,:));
h_ref   = r_n(3,:);

lat_mech = rad2deg(r_n_cal(1,:));
lon_mech = rad2deg(r_n_cal(2,:));
h_mech   = r_n_cal(3,:);

lat_ekf = rad2deg(r_n_ekf(1,:));
lon_ekf = rad2deg(r_n_ekf(2,:));
h_ekf   = r_n_ekf(3,:);

lat_fb = rad2deg(r_fb_ekf(1,:));
lon_fb = rad2deg(r_fb_ekf(2,:));
h_fb   = r_fb_ekf(3,:);

% Errors
lat_err_mech = lat_mech - lat_ref;
lon_err_mech = lon_mech - lon_ref;
h_err_mech   = h_mech   - h_ref;

lat_err_ekf = lat_ekf - lat_ref;
lon_err_ekf = lon_ekf - lon_ref;
h_err_ekf   = h_ekf   - h_ref;

lat_err_fb = lat_fb - lat_ref;
lon_err_fb = lon_fb - lon_ref;
h_err_fb   = h_fb   - h_ref;

figure('Name','Position - Comparison & Error','Color','w');

% --- Latitude ---
subplot(3,2,1);
plot(t, lat_mech, 'b',  'LineWidth', 1.2); hold on;
plot(t, lat_ekf,  'g--','LineWidth', 1.2);
plot(t, lat_fb,   'c-.','LineWidth', 1.2);
plot(t, lat_ref,  'r:', 'LineWidth', 1.2);
grid on;
legend('INS','EKF','FB EKF','Ref','Location','best');
title('Latitude');
ylabel('Latitude [deg]');
xlabel('Time [s]');

subplot(3,2,2);
plot(t, lat_err_mech, 'k', 'LineWidth', 1.2); hold on;
plot(t, lat_err_ekf,  'm', 'LineWidth', 1.2);
plot(t, lat_err_fb,   'c--','LineWidth', 1.2);
grid on;
legend('INS Error','EKF Error','FB Error','Location','best');
title('Latitude Error');
ylabel('Error [deg]');
xlabel('Time [s]');
ylim(autoErrorYLim([lat_err_mech, lat_err_ekf, lat_err_fb]));

% --- Longitude ---
subplot(3,2,3);
plot(t, lon_mech, 'b',  'LineWidth', 1.2); hold on;
plot(t, lon_ekf,  'g--','LineWidth', 1.2);
plot(t, lon_fb,   'c-.','LineWidth', 1.2);
plot(t, lon_ref,  'r:', 'LineWidth', 1.2);
grid on;
legend('INS','EKF','FB EKF','Ref','Location','best');
title('Longitude');
ylabel('Longitude [deg]');
xlabel('Time [s]');

subplot(3,2,4);
plot(t, lon_err_mech, 'k', 'LineWidth', 1.2); hold on;
plot(t, lon_err_ekf,  'm', 'LineWidth', 1.2);
plot(t, lon_err_fb,   'c--','LineWidth', 1.2);
grid on;
legend('INS Error','EKF Error','FB Error','Location','best');
title('Longitude Error');
ylabel('Error [deg]');
xlabel('Time [s]');
ylim(autoErrorYLim([lon_err_mech, lon_err_ekf, lon_err_fb]));

% --- Height ---
subplot(3,2,5);
plot(t, h_mech, 'b',  'LineWidth', 1.2); hold on;
plot(t, h_ekf,  'g--','LineWidth', 1.2);
plot(t, h_fb,   'c-.','LineWidth', 1.2);
plot(t, h_ref,  'r:', 'LineWidth', 1.2);
grid on;
legend('INS','EKF','FB EKF','Ref','Location','best');
title('Height');
ylabel('Height [m]');
xlabel('Time [s]');

subplot(3,2,6);
plot(t, h_err_mech, 'k', 'LineWidth', 1.2); hold on;
plot(t, h_err_ekf,  'm', 'LineWidth', 1.2);
plot(t, h_err_fb,   'c--','LineWidth', 1.2);
grid on;
legend('INS Error','EKF Error','FB Error','Location','best');
title('Height Error');
ylabel('Error [m]');
xlabel('Time [s]');
ylim(autoErrorYLim([h_err_mech, h_err_ekf, h_err_fb]));

%% ===================== EULER ANGLES (Roll, Pitch, Yaw) =============

% Wrapped angle errors in radians (mech, ekf, fb)
roll_err_mech  = angleError(roll_cal,  roll);
pitch_err_mech = angleError(pitch_cal, pitch);
yaw_err_mech   = angleError(yaw_cal,   yaw);

roll_err_ekf  = angleError(roll_ekf,  roll);
pitch_err_ekf = angleError(pitch_ekf, pitch);
yaw_err_ekf   = angleError(yaw_ekf,   yaw);

roll_err_fb  = angleError(roll_fb,  roll);
pitch_err_fb = angleError(pitch_fb, pitch);
yaw_err_fb   = angleError(yaw_fb,   yaw);

% Convert all to degrees
roll_ref_deg  = rad2deg(roll);
pitch_ref_deg = rad2deg(pitch);
yaw_ref_deg   = rad2deg(yaw);

roll_mech_deg  = rad2deg(roll_cal);
pitch_mech_deg = rad2deg(pitch_cal);
yaw_mech_deg   = rad2deg(yaw_cal);

roll_ekf_deg  = rad2deg(roll_ekf);
pitch_ekf_deg = rad2deg(pitch_ekf);
yaw_ekf_deg   = rad2deg(yaw_ekf);

roll_fb_deg  = rad2deg(roll_fb);
pitch_fb_deg = rad2deg(pitch_fb);
yaw_fb_deg   = rad2deg(yaw_fb);

roll_err_mech_deg  = rad2deg(roll_err_mech);
pitch_err_mech_deg = rad2deg(pitch_err_mech);
yaw_err_mech_deg   = rad2deg(yaw_err_mech);

roll_err_ekf_deg  = rad2deg(roll_err_ekf);
pitch_err_ekf_deg = rad2deg(pitch_err_ekf);
yaw_err_ekf_deg   = rad2deg(yaw_err_ekf);

roll_err_fb_deg  = rad2deg(roll_err_fb);
pitch_err_fb_deg = rad2deg(pitch_err_fb);
yaw_err_fb_deg   = rad2deg(yaw_err_fb);

figure('Name','Euler Angles - Comparison & Error','Color','w');

% --- Roll ---
subplot(3,2,1);
plot(t, roll_mech_deg, 'b',  'LineWidth', 1.2); hold on;
plot(t, roll_ekf_deg,  'g--','LineWidth', 1.2);
plot(t, roll_fb_deg,   'c-.','LineWidth', 1.2);
plot(t, roll_ref_deg,  'r:', 'LineWidth', 1.2);
grid on;
legend('INS','EKF','FB EKF','Ref','Location','best');
title('Roll');
ylabel('[deg]');
xlabel('Time [s]');

subplot(3,2,2);
plot(t, roll_err_mech_deg, 'k', 'LineWidth', 1.2); hold on;
plot(t, roll_err_ekf_deg,  'm', 'LineWidth', 1.2);
plot(t, roll_err_fb_deg,   'c--','LineWidth', 1.2);
grid on;
legend('INS Error','EKF Error','FB Error','Location','best');
title('Roll Error');
ylabel('Error [deg]');
xlabel('Time [s]');
ylim(autoErrorYLim([roll_err_mech_deg, roll_err_ekf_deg, roll_err_fb_deg]));

% --- Pitch ---
subplot(3,2,3);
plot(t, pitch_mech_deg, 'b',  'LineWidth', 1.2); hold on;
plot(t, pitch_ekf_deg,  'g--','LineWidth', 1.2);
plot(t, pitch_fb_deg,   'c-.','LineWidth', 1.2);
plot(t, pitch_ref_deg,  'r:', 'LineWidth', 1.2);
grid on;
legend('INS','EKF','FB EKF','Ref','Location','best');
title('Pitch');
ylabel('[deg]');
xlabel('Time [s]');

subplot(3,2,4);
plot(t, pitch_err_mech_deg, 'k', 'LineWidth', 1.2); hold on;
plot(t, pitch_err_ekf_deg,  'm', 'LineWidth', 1.2);
plot(t, pitch_err_fb_deg,   'c--','LineWidth', 1.2);
grid on;
legend('INS Error','EKF Error','FB Error','Location','best');
title('Pitch Error');
ylabel('Error [deg]');
xlabel('Time [s]');
ylim(autoErrorYLim([pitch_err_mech_deg, pitch_err_ekf_deg, pitch_err_fb_deg]));

% --- Yaw ---
subplot(3,2,5);
plot(t, yaw_mech_deg, 'b',  'LineWidth', 1.2); hold on;
plot(t, yaw_ekf_deg,  'g--','LineWidth', 1.2);
plot(t, yaw_fb_deg,   'c-.','LineWidth', 1.2);
plot(t, yaw_ref_deg,  'r:', 'LineWidth', 1.2);
grid on;
legend('INS','EKF','FB EKF','Ref','Location','best');
title('Yaw');
ylabel('[deg]');
xlabel('Time [s]');

subplot(3,2,6);
plot(t, yaw_err_mech_deg, 'k', 'LineWidth', 1.2); hold on;
plot(t, yaw_err_ekf_deg,  'm', 'LineWidth', 1.2);
plot(t, yaw_err_fb_deg,   'c--','LineWidth', 1.2);
grid on;
legend('INS Error','EKF Error','FB Error','Location','best');
title('Yaw Error');
ylabel('Error [deg]');
xlabel('Time [s]');
ylim(autoErrorYLim([yaw_err_mech_deg, yaw_err_ekf_deg, yaw_err_fb_deg]));

%% ===================== RMSE SUMMARY TABLE ==========================
% Compute RMSEs for all quantities

% Velocity [m/s]
rmse_vN_mech = computeRMSE(vel_err_mech(1,:));
rmse_vN_ekf  = computeRMSE(vel_err_ekf(1,:));
rmse_vN_fb   = computeRMSE(vel_err_fb(1,:));

rmse_vE_mech = computeRMSE(vel_err_mech(2,:));
rmse_vE_ekf  = computeRMSE(vel_err_ekf(2,:));
rmse_vE_fb   = computeRMSE(vel_err_fb(2,:));

rmse_vD_mech = computeRMSE(vel_err_mech(3,:));
rmse_vD_ekf  = computeRMSE(vel_err_ekf(3,:));
rmse_vD_fb   = computeRMSE(vel_err_fb(3,:));

% Position (lat,lon [deg], h [m])
rmse_lat_mech = computeRMSE(lat_err_mech);
rmse_lat_ekf  = computeRMSE(lat_err_ekf);
rmse_lat_fb   = computeRMSE(lat_err_fb);

rmse_lon_mech = computeRMSE(lon_err_mech);
rmse_lon_ekf  = computeRMSE(lon_err_ekf);
rmse_lon_fb   = computeRMSE(lon_err_fb);

rmse_h_mech = computeRMSE(h_err_mech);
rmse_h_ekf  = computeRMSE(h_err_ekf);
rmse_h_fb   = computeRMSE(h_err_fb);

% Attitude [deg]
rmse_roll_mech  = computeRMSE(roll_err_mech_deg);
rmse_roll_ekf   = computeRMSE(roll_err_ekf_deg);
rmse_roll_fb    = computeRMSE(roll_err_fb_deg);

rmse_pitch_mech = computeRMSE(pitch_err_mech_deg);
rmse_pitch_ekf  = computeRMSE(pitch_err_ekf_deg);
rmse_pitch_fb   = computeRMSE(pitch_err_fb_deg);

rmse_yaw_mech   = computeRMSE(yaw_err_mech_deg);
rmse_yaw_ekf    = computeRMSE(yaw_err_ekf_deg);
rmse_yaw_fb     = computeRMSE(yaw_err_fb_deg);

% Build table data
signalNames = { ...
    'v_N', 'v_E', 'v_D', ...
    'lat', 'lon', 'h', ...
    'roll', 'pitch', 'yaw'}';

units = { ...
    'm/s','m/s','m/s', ...
    'deg','deg','m', ...
    'deg','deg','deg'}';

rmseINS = [ ...
    rmse_vN_mech; rmse_vE_mech; rmse_vD_mech; ...
    rmse_lat_mech; rmse_lon_mech; rmse_h_mech; ...
    rmse_roll_mech; rmse_pitch_mech; rmse_yaw_mech];

rmseEKF = [ ...
    rmse_vN_ekf; rmse_vE_ekf; rmse_vD_ekf; ...
    rmse_lat_ekf; rmse_lon_ekf; rmse_h_ekf; ...
    rmse_roll_ekf; rmse_pitch_ekf; rmse_yaw_ekf];

rmseFB = [ ...
    rmse_vN_fb; rmse_vE_fb; rmse_vD_fb; ...
    rmse_lat_fb; rmse_lon_fb; rmse_h_fb; ...
    rmse_roll_fb; rmse_pitch_fb; rmse_yaw_fb];

% Create figure with uitable
figure('Name','RMSE Summary','Color','w');
data = [units, num2cell(rmseINS), num2cell(rmseEKF), num2cell(rmseFB)];
colNames = {'Unit','RMSE INS','RMSE EKF','RMSE FB'};
uitable('Data', data, ...
        'ColumnName', colNames, ...
        'RowName', signalNames, ...
        'Units','normalized', ...
        'Position',[0 0 1 1]);

end

%% ============== HELPER FUNCTIONS =====================

function err = angleError(a, b)
%ANGLEERROR  Wrapped angular difference: err = a - b in [-pi, pi]
    err = atan2(sin(a - b), cos(a - b));
end

function L = autoErrorYLim(y)
%AUTOERRORYLIM  Simple automatic y-limits for error plots.
    y = y(:);
    y = y(isfinite(y));  % remove NaN / Inf

    if isempty(y)
        L = [-1 1];
        return;
    end

    n = numel(y);
    if n < 2
        m = max(1, abs(y(1)));
        L = 1.1 * m * [-1 1];
        return;
    end

    % Use middle 90% of data (ignore first/last 5%)
    i1 = max(1, floor(0.05*n) + 1);
    i2 = max(i1 + 1, ceil(0.95*n));
    s  = y(i1:i2);

    s  = s(isfinite(s));
    if isempty(s)
        s = y;
    end

    lo = min(s);
    hi = max(s);

    if lo == hi
        span = max(1e-6, 0.1*max(1, abs(lo)));
        L = [lo - span, hi + span];
    else
        pad = 0.1 * (hi - lo);
        L = [lo - pad, hi + pad];
    end
end

function val = computeRMSE(e)
%COMPUTERMSE  Root-mean-square of a vector (ignores NaN/Inf)
    e = e(:);
    e = e(isfinite(e));
    if isempty(e)
        val = NaN;
    else
        val = sqrt(mean(e.^2));
    end
end
