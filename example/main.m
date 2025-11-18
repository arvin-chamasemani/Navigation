clc
clear
close all
% programmer: Arvin chamasrmani

% ensure script runs from its own folder
cd(fileparts(mfilename('fullpath')))

% loading data
load("../data/path2sensor_output.mat");

% main class
[r_cal, v_cal, q_cal, roll_cal, pitch_cal, yaw_cal] = ...
    Mechanization.runBatch(t, r_n, v_n, roll, pitch, yaw, f_b_IMU, w_b_ib_IMU);

% plotting
plotting(t, r_n, v_n, roll, pitch, yaw, ...
                r_cal, v_cal, roll_cal, pitch_cal, yaw_cal)

