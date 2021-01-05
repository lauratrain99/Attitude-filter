%% ATTITUDE COMPUTATION
% Author: Laura Train
% Date  21/12 .
%% 
% Kalman filter to estimate orientation simulating an IMU containing gyros
% and accelerometers. Simulate a static IMU with
% "Desarrollo de un Sistema Inercial de Referencia de Actitud basado en un
% Estimador Óptimo No Lineal"

%% Use NaveGo functions
matlabrc

addpath ../
addpath ../../ins/
addpath ../../simulation/
addpath ../../conversions/
addpath ../../kalman/
addpath ../../performance_analysis/

%% GENERATE SYNTHETIC DATA

imu1.ini_align_err = deg2rad([0.5, 0.5, 0.5]);
imu1.g_std = [0.01, 0.01, 0.01];
imu1.a_std = [0.01, 0.01, 0.01];
imu1.ini_align = deg2rad([0, 0, 0]);
imu1.gb_dyn = [0.001, 0.001, 0.001]; 

N = 10000;
imu1.t = linspace(0,10,N);
wbx = zeros(N,1) + imu1.g_std(1)*randn(N,1);
wby = zeros(N,1) + imu1.g_std(2)*randn(N,1); 
wbz = zeros(N,1) + imu1.g_std(3)*randn(N,1); 

imu1.wb = [wbx, wby, wbz];

abx = zeros(N,1) + imu1.a_std(1)*randn(N,1);
aby = zeros(N,1) + imu1.a_std(2)*randn(N,1); 
abz = -9.81*ones(N,1) + imu1.a_std(3)*randn(N,1);

imu1.fb = [abx, aby, abz];

% mnx = 0.22*ones(length(imu1.t),1);
% mny = zeros(length(imu1.t),1);
% mnz = 0.17*ones(length(imu1.t),1);
% 
% imu1.mn = [mnx, mny, mnz];

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

% figure(3)
% plot(imu1.t,imu1.mn(:,1), 'r', imu1.t, imu1.mn(:,2), 'b', imu1.t, imu1.mn(:,3), 'g')
% xlabel('Time [s]')
% grid on
% title('IMU1 Local magnetic field [Gauss]')
% legend('mx','my','mz')
% saveas(figure(3),'IMU1_magnetometer_raw.jpg')

%% SENSOR FUSION AND FILTER

[nav1, kf1] = imu_filter(imu1);



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