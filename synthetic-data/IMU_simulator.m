function [imu] = IMU_simulator(imu)
   
    N = length(imu.t);
    [~, euler] = attitude_computer(imu);

    roll = euler(:,1);
    pitch = euler(:,2);
    yaw = euler(:,3);
    
    for i=1:N
        imu.fb(i,1:3)=euler2dcm([roll(i), pitch(i), yaw(i)])*[0;0;-9.8];
    end


    for i=1:N
        imu.mb(i,1:3)=euler2dcm([roll(i), pitch(i), yaw(i)])*[0.22;0;0.17];
    end
end


% data = [t', wb, fb, mb];
% fileID = fopen('xAxisRotation.txt','w');
% fprintf(fileID,'%6.2f %6.2f %12.8f %6.2f %12.8f %6.2f %12.8f %6.2f %6.2f %12.8f\n',data');
% fclose(fileID);

 
% mag = [0.22; 0; 1.17];
% for i=1:length(t)
% 
%     cmTerr(1:3,i)=euler2dcm([roll(i), pitch(i), yaw(i)])'*[0.22; 0; 1.17];
% 
% 
% end