%% ATTITUDE COMPUTATION
% Author: Laura Train
% Date 04/01
%% 
% Kalman filter to estimate orientation simulating an IMU containing three-
% axis gyro, accelerometer and magnetometer
% Constant angular velocity in x axis

%% Use NaveGo functions
matlabrc

addpath ../
addpath ../../ins/
addpath ../../simulation/
addpath ../../conversions/
addpath ../../kalman/

%% GENERATE SYNTHETIC DATA

imu1.ini_align_err = deg2rad([0.5, 0.5, 0.5]);
imu1.g_std = [0.05, 0.05, 0.05];
imu1.a_std = [0.01, 0.01, 0.01];
imu1.m_std = [0.002, 0.002, 0.002];
imu1.ini_align = deg2rad([1, 2, 3]);
imu1.gb_dyn = [0.001, 0.001, 0.001]; 

imu1.t = 0:1/100:10;
N = length(imu1.t);

imu1.wb(:,1) = 4*pi/10*ones(N,1);
imu1.wb(:,2) = zeros(N, 1);
imu1.wb(:,3) = zeros(N, 1);

[imu1] = IMU_simulator(imu1);

imu1.wb(:,1) = imu1.wb(:,1) + imu1.g_std(1)*randn(N,1);
imu1.wb(:,2) = imu1.wb(:,2) + imu1.g_std(2)*randn(N,1);
imu1.wb(:,3) = imu1.wb(:,3) + imu1.g_std(3)*randn(N,1);
imu1.fb(:,1) = imu1.fb(:,1) + imu1.a_std(1)*randn(N,1);
imu1.fb(:,2) = imu1.fb(:,2) + imu1.a_std(2)*randn(N,1);
imu1.fb(:,3) = imu1.fb(:,3) + imu1.a_std(3)*randn(N,1);
imu1.mb(:,1) = imu1.mb(:,1) + imu1.m_std(1)*randn(N,1);
imu1.mb(:,2) = imu1.mb(:,2) + imu1.m_std(2)*randn(N,1);
imu1.mb(:,3) = imu1.mb(:,3) + imu1.m_std(3)*randn(N,1);

%% PLOT SYNTHETIC DATA: GYRO + ACC + MAG

figure(1)
plot(imu1.t,imu1.fb(:,1), 'r', imu1.t, imu1.fb(:,2), 'b', imu1.t, imu1.fb(:,3), 'g')
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')
grid on
title('IMU1 Acceleration raw measurements body frame [m/s^2]')
legend('ax','ay','az')
saveas(figure(1),'IMU1_accelerations_raw.jpg')

figure(2)
plot(imu1.t,imu1.wb(:,1), 'r', imu1.t, imu1.wb(:,2), 'b', imu1.t, imu1.wb(:,3), 'g')
xlabel('Time [s]')
ylabel('Angular velocity [rad/s]')
grid on
title('IMU1 Angular velocity raw measurements body frame [rad/s]')
legend('wx','wy','wz')
saveas(figure(2),'IMU1_angularvelocities_raw.jpg')

figure(3)
plot(imu1.t,imu1.wb(:,1), 'r', imu1.t, imu1.wb(:,2), 'b', imu1.t, imu1.wb(:,3), 'g')
xlabel('Time [s]')
ylabel('Magnetic field [Gauss]')
grid on
title('IMU1 Local magnetic field raw measurements body frame [Gauss]')
legend('wx','wy','wz')
saveas(figure(2),'IMU1_magneticfield_raw.jpg')

%% SENSOR FUSION AND FILTER

[nav1, kf1] = imu_filter_magnetometer(imu1);

%% COMPARISON ATTITUDE COMPUTER

[quat, euler] = attitude_computer(imu1);

for i = 1:length(nav1.t)
    nav1.quat_error_norm(i,1) = norm(nav1.deltaxp(1:3,:));
    nav1.dyn_bias_norm(i,1) = norm(nav1.deltaxp(4:6,:));
end


figure(3)
plot(nav1.t, nav1.quat_error_norm, 'r', nav1.t, nav1.dyn_bias_norm, 'b')
xlabel('Time [s]')
grid minor
legend('\deltaq','\delta\zeta')
title('Errors')

figure(4)
plot(imu1.t, quat(:,1), 'r', imu1.t, quat(:,2), 'c', imu1.t, quat(:,3), 'g', imu1.t, quat(:,4), 'k', ...
     nav1.t, nav1.qua(:,1), 'or', nav1.t, nav1.qua(:,2), 'oc', nav1.t, nav1.qua(:,3), 'og', nav1.t, nav1.qua(:,4), 'ok')
xlabel('Time [s]')
legend('q1', 'q2', 'q3', 'q4','q1 Kalman','q2 Kalman','q3 Kalman', 'q4 Kalman')
grid minor
title('Attitude computer vs Kalman filter. Quaternions')
legend('location','southeast')

figure(5)
plot(imu1.t, rad2deg(euler(:,1)), 'r', imu1.t, rad2deg(euler(:,2)), 'b', imu1.t, rad2deg(euler(:,3)), 'g', ...
     nav1.t, nav1.roll, 'or', nav1.t, nav1.pitch, 'ob', nav1.t, nav1.yaw, 'og')
legend('roll', 'pitch', 'yaw', 'roll Kalman', 'pitch Kalman', 'yaw Kalman')
xlabel('Time [s]')
grid minor
title('Attitude computer vs Kalman filter. Euler angles [deg]')

