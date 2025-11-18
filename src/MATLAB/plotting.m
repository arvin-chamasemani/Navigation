function plotting(t, r_n, v_n, roll, pitch, yaw, ...
                  r_n_cal, v_n_cal, roll_cal, pitch_cal, yaw_cal)
%PLOT_INS  Compare reference vs. calculated INS states (velocity, position, attitude) and plot their errors.

%% ===================== VELOCITY (N,E,D) ============================
% Errors
vel_err = v_n_cal - v_n;

figure('Name','Velocity (NED) - Comparison & Error','Color','w');

% --- North velocity ---
subplot(3,2,1);
plot(t, v_n_cal(1,:), 'LineWidth', 1.2); hold on;
plot(t, v_n(1,:), '--r', 'LineWidth', 1.0);
grid on;
legend('Calc','Ref','Location','best');
title('North Velocity');
ylabel('v_N [m/s]');
xlabel('Time [s]');

subplot(3,2,2);
plot(t, vel_err(1,:), 'k', 'LineWidth', 1.2);
grid on;
title('North Velocity Error');
ylabel('Error [m/s]');
xlabel('Time [s]');
ylim(autoErrorYLim(vel_err(1,:)));

% --- East velocity ---
subplot(3,2,3);
plot(t, v_n_cal(2,:), 'LineWidth', 1.2); hold on;
plot(t, v_n(2,:), '--r', 'LineWidth', 1.0);
grid on;
legend('Calc','Ref','Location','best');
title('East Velocity');
ylabel('v_E [m/s]');
xlabel('Time [s]');

subplot(3,2,4);
plot(t, vel_err(2,:), 'k', 'LineWidth', 1.2);
grid on;
title('East Velocity Error');
ylabel('Error [m/s]');
xlabel('Time [s]');
ylim(autoErrorYLim(vel_err(2,:)));

% --- Down velocity ---
subplot(3,2,5);
plot(t, v_n_cal(3,:), 'LineWidth', 1.2); hold on;
plot(t, v_n(3,:), '--r', 'LineWidth', 1.0);
grid on;
legend('Calc','Ref','Location','best');
title('Down Velocity');
ylabel('v_D [m/s]');
xlabel('Time [s]');

subplot(3,2,6);
plot(t, vel_err(3,:), 'k', 'LineWidth', 1.2);
grid on;
title('Down Velocity Error');
ylabel('Error [m/s]');
xlabel('Time [s]');
ylim(autoErrorYLim(vel_err(3,:)));

%% ===================== POSITION (Lat, Lon, H) ======================

% Convert lat/lon to degrees
lat_ref = rad2deg(r_n(1,:));
lon_ref = rad2deg(r_n(2,:));
h_ref   = r_n(3,:);

lat_cal = rad2deg(r_n_cal(1,:));
lon_cal = rad2deg(r_n_cal(2,:));
h_cal   = r_n_cal(3,:);

% Errors
lat_err = lat_cal - lat_ref;
lon_err = lon_cal - lon_ref;
h_err   = h_cal   - h_ref;

figure('Name','Position - Comparison & Error','Color','w');

% --- Latitude ---
subplot(3,2,1);
plot(t, lat_cal, 'LineWidth', 1.2); hold on;
plot(t, lat_ref, '--r', 'LineWidth', 1.0);
grid on;
legend('Calc','Ref','Location','best');
title('Latitude');
ylabel('Latitude [deg]');
xlabel('Time [s]');

subplot(3,2,2);
plot(t, lat_err, 'k', 'LineWidth', 1.2);
grid on;
title('Latitude Error');
ylabel('Error [deg]');
xlabel('Time [s]');
ylim(autoErrorYLim(lat_err));

% --- Longitude ---
subplot(3,2,3);
plot(t, lon_cal, 'LineWidth', 1.2); hold on;
plot(t, lon_ref, '--r', 'LineWidth', 1.0);
grid on;
legend('Calc','Ref','Location','best');
title('Longitude');
ylabel('Longitude [deg]');
xlabel('Time [s]');

subplot(3,2,4);
plot(t, lon_err, 'k', 'LineWidth', 1.2);
grid on;
title('Longitude Error');
ylabel('Error [deg]');
xlabel('Time [s]');
ylim(autoErrorYLim(lon_err));

% --- Height ---
subplot(3,2,5);
plot(t, h_cal, 'LineWidth', 1.2); hold on;
plot(t, h_ref, '--r', 'LineWidth', 1.0);
grid on;
legend('Calc','Ref','Location','best');
title('Height');
ylabel('Height [m]');
xlabel('Time [s]');

subplot(3,2,6);
plot(t, h_err, 'k', 'LineWidth', 1.2);
grid on;
title('Height Error');
ylabel('Error [m]');
xlabel('Time [s]');
ylim(autoErrorYLim(h_err));

%% ===================== EULER ANGLES (Roll, Pitch, Yaw) =============

% Wrap angle errors in radians
roll_err  = angleError(roll_cal,  roll);
pitch_err = angleError(pitch_cal, pitch);
yaw_err   = angleError(yaw_cal,   yaw);

% Convert all to degrees using rad2deg
roll_ref_deg  = rad2deg(roll);
pitch_ref_deg = rad2deg(pitch);
yaw_ref_deg   = rad2deg(yaw);

roll_cal_deg  = rad2deg(roll_cal);
pitch_cal_deg = rad2deg(pitch_cal);
yaw_cal_deg   = rad2deg(yaw_cal);

roll_err_deg  = rad2deg(roll_err);
pitch_err_deg = rad2deg(pitch_err);
yaw_err_deg   = rad2deg(yaw_err);

figure('Name','Euler Angles - Comparison & Error','Color','w');

% --- Roll ---
subplot(3,2,1);
plot(t, roll_cal_deg, 'LineWidth', 1.2); hold on;
plot(t, roll_ref_deg, '--r', 'LineWidth', 1.0);
grid on;
legend('Calc','Ref','Location','best');
title('Roll');
ylabel('[deg]');
xlabel('Time [s]');

subplot(3,2,2);
plot(t, roll_err_deg, 'k', 'LineWidth', 1.2);
grid on;
title('Roll Error');
ylabel('Error [deg]');
xlabel('Time [s]');
ylim(autoErrorYLim(roll_err_deg));

% --- Pitch ---
subplot(3,2,3);
plot(t, pitch_cal_deg, 'LineWidth', 1.2); hold on;
plot(t, pitch_ref_deg, '--r', 'LineWidth', 1.0);
grid on;
legend('Calc','Ref','Location','best');
title('Pitch');
ylabel('[deg]');
xlabel('Time [s]');

subplot(3,2,4);
plot(t, pitch_err_deg, 'k', 'LineWidth', 1.2);
grid on;
title('Pitch Error');
ylabel('Error [deg]');
xlabel('Time [s]');
ylim(autoErrorYLim(pitch_err_deg));

% --- Yaw ---
subplot(3,2,5);
plot(t, yaw_cal_deg, 'LineWidth', 1.2); hold on;
plot(t, yaw_ref_deg, '--r', 'LineWidth', 1.0);
grid on;
legend('Calc','Ref','Location','best');
title('Yaw');
ylabel('[deg]');
xlabel('Time [s]');

subplot(3,2,6);
plot(t, yaw_err_deg, 'k', 'LineWidth', 1.2);
grid on;
title('Yaw Error');
ylabel('Error [deg]');
xlabel('Time [s]');
ylim(autoErrorYLim(yaw_err_deg));

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
