%% ATTITUDE COMPUTATION
% Author: Laura Train
% Date 15/11
% The goal is to put into practice the theory explained in the book 
% "Desarrollo de un Sistema Inercial de Referencia de Actitud basado en un
% Estimador Óptimo No Lineal", specifically section two: 
% Attitude Computation, Numerical Integration

%% Use NaveGo functions
matlabrc

addpath ../../ins/
addpath ../../ins-gnss/
addpath ../../simulation/
addpath ../../conversions/
addpath ../../performance_analysis/

%% CASE 2: wx, wy follow trigonometric functions and wz is constant
t = linspace(0,10,361)';
std_w = [0.1, 0.1, 0.05];

wbx = sin(t) + std_w(1)*randn(length(t),1);
wby = cos(t) + std_w(2)*randn(length(t),1);
wbz = 1.5*ones(length(t),1) + std_w(3)*randn(length(t),1); 

wb = [wbx, wby, wbz];

figure(1)
plot(t, wbx, 'r', t, wby, 'b', t, wbz, 'g')
grid minor
xlabel('time [s]')
ylabel('Angular velocities[rad/s]')
legend('wx','wy','wz')
legend('location','northeast')
title('Gyroscope synthetic data [body frame]')


wnorm_initial = sqrt(wbx(1) + wby(1) + wbz(1));

if wnorm_initial == 0
    q_initial = [0, 0, 0, 1];
else
    dt = t(2) - t(1);
    n1 = wbx(1)/wnorm_initial;
    n2 = wby(2)/wnorm_initial;
    n3 = wbz(3)/wnorm_initial;
    
    q_initial(1) = n1*sin(0.5*wnorm_initial*dt);
    q_initial(2) = n2*sin(0.5*wnorm_initial*dt);
    q_initial(3) = n3*sin(0.5*wnorm_initial*dt);
    q_initial(4) = cos(0.5*wnorm_initial*dt);
    % no need of normalizing. Already normalized
end

qua = q_initial;
% [a11, a21, a31, a21, a22, a23, a31, a32, a33]
DCMbn_n = reshape(qua2dcm(qua),1, 9);               
euler   = rad2deg(qua2euler(qua));

for i = 2:length(t)
    dt = t(i) - t(i-1);
    [qua(i,:), DCMbn_n(i,:), euler(i,:)] = my_quat_update(wb(i,:), qua(i-1,:)', dt);
end

roll = rad2deg(unwrap(euler(:,3)));
pitch = rad2deg(unwrap(euler(:,2)));
yaw = rad2deg(unwrap(euler(:,1)));

figure(2)
plot(t, roll, 'r', t, pitch, 'b', t, yaw, 'g')
xlabel('Time [s]')
legend('roll','pitch','yaw')
title('Euler angles [º]')
grid minor
legend('location','southwest')

q1 = qua(:,1);
q2 = qua(:,2);
q3 = qua(:,3);
q4 = qua(:,4);

figure(3)
plot(t, q1, 'r', t, q2, 'b', t, q3, 'g', t, q4, 'k')
xlabel('Time [s]')
legend('q1','q2','q3','q4')
title('Quaternions')
grid minor
legend('location','southwest')

%% CASE 3:angular velocity has a fixed direction
matlabrc
freq_gyro = 100;
t = linspace(0, 5, freq_gyro)';

freq_ex = 1;
a = 1;
Omega = 2*pi*freq_ex;

% Exact solution
omegai1 = sin(a)*sin(Omega*t);
omegai2 = zeros(length(t),1);
omegai3 = zeros(length(t),1); 

omegai = [omegai1, omegai2, omegai3];

q_initial = [0, 0, 0, 1];

qua_real = q_initial;

for i = 2:length(t)
    w = 1/2 * omegai(i,:);
    wnorm = norm(w);

    if wnorm == 0

        cos_const = 1;
        sin_const = 1;
    else

        cos_const=cos(wnorm);
        sin_const=sin(wnorm)/wnorm;
    end

    W=[0  -w(1) -w(2) -w(3);
     w(1)    0  -w(3)  w(2);
     w(2)  w(3)    0  -w(1);
     w(3) -w(2)  w(1)     0];
    I = eye(4);
    qua_real(i,:) = (I*cos_const + W*sin_const) * q_initial';

end

q1_real = qua_real(:,1);
q2_real = qua_real(:,2);
q3_real = qua_real(:,3);
q4_real = qua_real(:,4);

figure(4)
plot(t, q1_real, 'r', t, q2_real, 'b', t, q3_real, 'g', t, q4_real, 'k')
xlabel('Time [s]')
grid minor
title('Quaternions for constant direction angular velocity')

% Numerical solution 100 Hz

omega1 = Omega*sin(a)*cos(Omega*t);
omega2 = zeros(length(t),1);
omega3 = zeros(length(t),1); 

omega = [omega1, omega2, omega3];
qua = q_initial;

for i = 2:length(t)
    dt = t(i) - t(i-1);
    [qua(i,:), ~, ~] = my_quat_update(omega(i,:), qua(i-1,:)', dt);
end

q1 = qua(:,1);
q2 = qua(:,2);
q3 = qua(:,3);
q4 = qua(:,4);

hold on
plot(t, q1, '--r', t, q2, '--b', t, q3, '--g', t, q4, '--k')
legend('q1 real','q2 real','q3 real','q4 real','q1 100 Hz','q2 100 Hz','q3 100 Hz','q4 100 Hz')
legend('location','best')

%% CASE 3: conical motion
matlabrc
freq_gyro = 100;
t = linspace(0, 5, freq_gyro)';

freq_ex = 1;
a = 1;
Omega = 2*pi*freq_ex;

omega1 = 2*Omega*sin(a/2)^2*ones(length(t),1);
omega2 = -Omega*sin(a)*sin(Omega*t);
omega3 = Omega*sin(a)*cos(Omega*t); 

omega = [omega1, omega2, omega3];

% Exact solution
q1_real = zeros(length(t),1);
q2_real = sin(a/2)*cos(Omega*t);
q3_real = sin(a/2)*sin(Omega*t);
q4_real = cos(a/2)*ones(length(t),1);


figure(5)
plot(t, q1_real, 'r', t, q2_real, 'b', t, q3_real, 'g', t, q4_real, 'k')
xlabel('Time [s]')
grid minor
title('Quaternions for pure cone motion')


% Numeric solution 100 Hz
q_initial = [0, sin(a/2), 0, cos(a/2)];

qua = q_initial;


for i = 2:length(t)
    dt = t(i) - t(i-1);
    [qua(i,:), ~, ~] = my_quat_update(omega(i,:), qua(i-1,:)', dt);
end


q1 = qua(:,1);
q2 = qua(:,2);
q3 = qua(:,3);
q4 = qua(:,4);

hold on
plot(t, q1, '--r', t, q2, '--b', t, q3, '--g', t, q4, '--k')
legend('q1 real','q2 real','q3 real','q4 real','q1 100 Hz','q2 100 Hz','q3 100 Hz','q4 100 Hz')
legend('location','best')