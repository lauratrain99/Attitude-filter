%% ATTITUDE COMPUTATION
% Author: Laura Train
% Date 27/12
%% 
% Kalman filter to estimate orientation simulating an IMU containing gyros
% and accelerometers 
% Constat angular velocity in x axis

%% Use NaveGo functions
matlabrc

addpath ../
addpath ../../ins/
addpath ../../simulation/
addpath ../../conversions/
addpath ../../kalman/
addpath ../../performance_analysis/

%% CONVERSION CONSTANTS

G =  9.80665;       % Gravity constant, m/s^2
G2MSS = G;          % g to m/s^2
MSS2G = (1/G);      % m/s^2 to g

D2R = (pi/180);     % degrees to radians
R2D = (180/pi);     % radians to degrees

KT2MS = 0.514444;   % knot to m/s
MS2KMH = 3.6;       % m/s to km/h
params = gyroparams;

%% REFERENCE DATA

load ref.mat

% ref.mat contains a reference data structure from which inertial
% sensors wil be simulated. It must contain the following fields:

%         t: Nx1 time vector (seconds).
%       lat: Nx1 latitude (radians).
%       lon: Nx1 longitude (radians).
%         h: Nx1 altitude (m).
%       vel: Nx3 NED velocities (m/s).
%      roll: Nx1 roll angles (radians).
%     pitch: Nx1 pitch angles (radians).
%       yaw: Nx1 yaw angle vector (radians).
%   DCMnb_m: Nx9 matrix with nav-to-body direct cosine matrices (DCM).
%            Each row of DCMnb_m contains the 9 elements of a particular DCMnb
%            matrix ordered as [a11 a21 a31 a12 a22 a32 a13 a23 a33].
%      freq: sampling frequency (Hz).

%% ADIS16405: IMU1 error profile

% IMU data structure:
%         t: Ix1 time vector (seconds).
%        fb: Ix3 accelerations vector in body frame XYZ (m/s^2).
%        wb: Ix3 turn rates vector in body frame XYZ (radians/s).
%       arw: 1x3 angle random walks (rad/s/root-Hz).
%      arrw: 1x3 angle rate random walks (rad/s^2/root-Hz).
%       vrw: 1x3 velocity random walks (m/s^2/root-Hz).
%      vrrw: 1x3 velocity rate random walks (m/s^3/root-Hz).
%     g_std: 1x3 gyros standard deviations (radians/s).
%     a_std: 1x3 accrs standard deviations (m/s^2).
%    gb_sta: 1x3 gyros static biases or turn-on biases (radians/s).
%    ab_sta: 1x3 accrs static biases or turn-on biases (m/s^2).
%    gb_dyn: 1x3 gyros dynamic biases or bias instabilities (radians/s).
%    ab_dyn: 1x3 accrs dynamic biases or bias instabilities (m/s^2).
%   gb_corr: 1x3 gyros correlation times (seconds).
%   ab_corr: 1x3 accrs correlation times (seconds).
%    gb_psd: 1x3 gyros dynamic biases PSD (rad/s/root-Hz).
%    ab_psd: 1x3 accrs dynamic biases PSD (m/s^2/root-Hz);
%      freq: 1x1 sampling frequency (Hz).
% ini_align: 1x3 initial attitude at t(1), [roll pitch yaw] (rad).
% ini_align_err: 1x3 initial attitude errors at t(1), [roll pitch yaw] (rad).

ADIS16405.arw      = 2   .* ones(1,3);     % Angle random walks [X Y Z] (deg/root-hour)
ADIS16405.arrw     = zeros(1,3);           % Angle rate random walks [X Y Z] (deg/root-hour/s)
ADIS16405.vrw      = 0.2 .* ones(1,3);     % Velocity random walks [X Y Z] (m/s/root-hour)
ADIS16405.vrrw     = zeros(1,3);           % Velocity rate random walks [X Y Z] (deg/root-hour/s)
ADIS16405.gb_sta   = 3   .* ones(1,3);     % Gyro static biases [X Y Z] (deg/s)
ADIS16405.ab_sta   = 50  .* ones(1,3);     % Acc static biases [X Y Z] (mg)
ADIS16405.gb_dyn   = 0.007 .* ones(1,3);   % Gyro dynamic biases [X Y Z] (deg/s)
ADIS16405.ab_dyn   = 0.2 .* ones(1,3);     % Acc dynamic biases [X Y Z] (mg)
ADIS16405.gb_corr  = 100 .* ones(1,3);     % Gyro correlation times [X Y Z] (seconds)
ADIS16405.ab_corr  = 100 .* ones(1,3);     % Acc correlation times [X Y Z] (seconds)
ADIS16405.freq     = ref.freq;             % IMU operation frequency [X Y Z] (Hz)
ADIS16405.m_psd    = 0.066 .* ones(1,3);   % Magnetometer noise density [X Y Z] (mgauss/root-Hz)

% ref time is used to simulate IMU sensors
ADIS16405.t = ref.t;                       % IMU time vector
dt = mean(diff(ADIS16405.t));              % IMU sampling interval

imu1 = imu_si_errors(ADIS16405, dt);       % IMU manufacturer error units to SI units.

imu1.ini_align_err = [3 3 10] .* D2R;                   % Initial attitude align errors for matrix P in Kalman filter, [roll pitch yaw] (radians)
imu1.ini_align = [0 0 0]; % Initial attitude align at t(1) (radians).

%% GENERATE SYNTHETIC DATA
% Generate N samples at a sampling rate of Fs with a sinusoidal frequency
% of Fc.

load('xAxisRotation.mat')

%%
imu1.t = xAxisRotation(:,1);
N = length(imu1.t);
imu1.wb(:,1) = xAxisRotation(:,2) + imu1.g_std(1)*randn(N,1);
imu1.wb(:,2) = xAxisRotation(:,3) + imu1.g_std(2)*randn(N,1);
imu1.wb(:,3) = xAxisRotation(:,4) + imu1.g_std(3)*randn(N,1);
imu1.fb(:,1) = xAxisRotation(:,5) + imu1.a_std(1)*randn(N,1);
imu1.fb(:,2) = xAxisRotation(:,6) + imu1.a_std(2)*randn(N,1);
imu1.fb(:,3) = xAxisRotation(:,7) + imu1.a_std(3)*randn(N,1);





figure(1)
plot(imu1.t,imu1.fb(:,1), 'r', imu1.t, imu1.fb(:,2), 'b', imu1.t, imu1.fb(:,3), 'g')
xlabel('Time [s]')
grid on
title('IMU1 Acceleration raw measurements body frame [m/s^2]')
legend('ax','ay','az')
saveas(figure(1),'IMU1_accelerations_raw.jpg')

figure(2)
plot(imu1.t,imu1.wb(:,1), 'r', imu1.t, imu1.wb(:,2), 'b', imu1.t, imu1.wb(:,3), 'g')
xlabel('Time [s]')
grid on
title('IMU1 Angular velocity measurements body frame [rad/s]')
legend('wx','wy','wz')
saveas(figure(2),'IMU1_angularvelocities_raw.jpg')

%% SENSOR FUSION AND FILTER

[nav1, kf1] = imu_filter(imu1);

figure(3)
plot(nav1.t, nav1.qua(:,1), 'r', nav1.t, nav1.qua(:,2), 'c', nav1.t, nav1.qua(:,3), 'g', nav1.t, nav1.qua(:,4), 'k')
xlabel('Time [s]')
legend('q1', 'q2', 'q3', 'q4')
grid minor
title('ADIS16405 (IMU1) quaternions')
legend('location','southeast')
saveas(figure(3),'IMU1_quaternions_EKF.jpg')


figure(4)
plot(nav1.t, nav1.roll, 'r', nav1.t, nav1.pitch, 'b', nav1.t, nav1.yaw, 'g')
legend('roll [deg]', 'pitch [deg]', 'yaw [deg]')
xlabel('Time [s]')
grid minor
title('ADIS16405 (IMU1) Euler angles')
saveas(figure(4),'IMU1_EulerAngles_EKF.jpg')

%% COMPARISON ATTITUDE COMPUTER
[quat, euler] = attitude_computer(imu1);


figure(5)
plot(imu1.t, quat(:,1), 'r', imu1.t, quat(:,2), 'c', imu1.t, quat(:,3), 'g', imu1.t, quat(:,4), 'k')
xlabel('Time [s]')
legend('q1', 'q2', 'q3', 'q4')
grid minor
title('ADIS16405 (IMU1) attitude computer quaternions')
legend('location','southeast')



figure(6)
plot(imu1.t, rad2deg(euler(:,1)), 'r', imu1.t, rad2deg(euler(:,2)), 'b', imu1.t, rad2deg(euler(:,3)), 'g')
legend('roll [deg]', 'pitch [deg]', 'yaw [deg]')
xlabel('Time [s]')
grid minor
title('ADIS16405 (IMU1) attitude computer Euler angles')

for i = 1:length(nav1.t)
    nav1.quat_error_norm(i,1) = norm(nav1.deltaxp(1:3,:));
    nav1.dyn_bias_norm(i,1) = norm(nav1.deltaxp(4:6,:));
end


figure(7)
plot(nav1.t, nav1.quat_error_norm, 'r', nav1.t, nav1.dyn_bias_norm, 'b')
xlabel('Time [s]')
grid minor
legend('\deltaq','\delta\zeta')
%xlim([-20, nav1.t(end) + 20])
title('Errors')


figure(8)
plot(imu1.t, quat(:,1), 'r', imu1.t, quat(:,2), 'c', imu1.t, quat(:,3), 'g', imu1.t, quat(:,4), 'k', ...
     nav1.t, nav1.qua(:,1), '--r', nav1.t, nav1.qua(:,2), '--c', nav1.t, nav1.qua(:,3), '--g', nav1.t, nav1.qua(:,4), '--k')
xlabel('Time [s]')
legend('q1', 'q2', 'q3', 'q4','q1 Kalman','q2 Kalman','q3 Kalman', 'q4 Kalman')
grid minor
title('Attitude computer vs Kalman filter. Quaternions')
legend('location','southeast')

figure(9)
plot(imu1.t, rad2deg(euler(:,1)), 'r', imu1.t, rad2deg(euler(:,2)), 'b', imu1.t, rad2deg(euler(:,3)), 'g', ...
     nav1.t, nav1.roll, '--r', nav1.t, nav1.pitch, '--b', nav1.t, nav1.yaw, '--g')
legend('roll', 'pitch', 'yaw', 'roll Kalman', 'pitch Kalman', 'yaw Kalman')
xlabel('Time [s]')
grid minor
title('Attitude computer vs Kalman filter. Euler angles [deg]')