% 比较两个惯导器件的输出
% IMU5210传输数据要5ms，ADIS16448传输数据要350us，IMU5210有时延5ms
% IMU5210大约有100ms的滤波延时

% 以各自的时间戳为x轴坐标，减去整数公共时间，让时间从0开始
t0 = floor(min([IMU5210_data(1,1),ADIS16448_data(1,1)]));
t_IMU5210 = IMU5210_data(:,1) - t0;
t_ADIS16448 = ADIS16448_data(:,1) - t0;

figure
subplot(3,2,1)
plot(t_IMU5210,IMU5210_data(:,2))
hold on
plot(t_ADIS16448,ADIS16448_data(:,2))

subplot(3,2,3)
plot(t_IMU5210,IMU5210_data(:,3))
hold on
plot(t_ADIS16448,ADIS16448_data(:,3))

subplot(3,2,5)
plot(t_IMU5210,IMU5210_data(:,4))
hold on
plot(t_ADIS16448,ADIS16448_data(:,4))

subplot(3,2,2)
plot(t_IMU5210,IMU5210_data(:,5))
hold on
plot(t_ADIS16448,ADIS16448_data(:,5))

subplot(3,2,4)
plot(t_IMU5210,IMU5210_data(:,6))
hold on
plot(t_ADIS16448,ADIS16448_data(:,6))

subplot(3,2,6)
plot(t_IMU5210,IMU5210_data(:,7))
hold on
plot(t_ADIS16448,ADIS16448_data(:,7))