function [quat, euler] = attitude_computer(imu)
    wb = imu.wb;
    t = imu.t;
    
    % initialization
    euler = [imu.ini_align(1), imu.ini_align(2), imu.ini_align(3)];
    quat = euler2qua(euler)';
    
    % [a11, a21, a31, a21, a22, a23, a31, a32, a33]
    for i = 2:length(t)
        dt = t(i) - t(i-1);
        [quat(i,:), ~, euler(i,:)] = my_quat_update(wb(i,:), quat(i-1,:)', dt);
    end
end

