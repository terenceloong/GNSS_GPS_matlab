% 解析IMU数据，将时间戳转化为GPS时间

%% 1.读文件
filename = '.\data\7_2\IMU5210_20190702.dat';
fileID = fopen(filename, 'r');
stream = fread(fileID, 'uint8=>uint8');
fclose(fileID);
device = strtok(filename, '_'); %设备名

%% 2.解析原始数据
n = length(stream); %总字节数
data = zeros(ceil(n/27),16); %解析的原始数据，每帧27字节
% [cnt, year,mon,day, hour,min,sec, TIM3, wx,wy,wz, fx,fy,fz, temp, PPS_error]

k = 1; %字节指针
m = 1; %数据存储指针
while 1
    if stream(k)==85 %帧头，0x55
        if k+26<=n %帧尾在总字节数以内
            if stream(k+26)==170 %帧尾，0xAA
                buff = stream(k+(0:26)); %提取一帧
                data(m,1) = double(buff(3)); %cnt
                data(m,2) = double(buff(4)); %year
                data(m,3) = double(buff(5)); %mon
                data(m,4) = double(buff(6)); %day
                data(m,5) = double(buff(7)); %hour
                data(m,6) = double(buff(8)); %min
                data(m,7) = double(buff(9)); %sec
                data(m,8)  = double(typecast(buff(10:11),'uint16')); %TIM3
                switch buff(2) %根据设备号进行数据转换
                    case 0
                        data(m,10) =  double(typecast(buff(12:13),'int16')) /32768*300;
                        data(m,9)  =  double(typecast(buff(14:15),'int16')) /32768*300;
                        data(m,11) = -double(typecast(buff(16:17),'int16')) /32768*300;
                        data(m,13) =  double(typecast(buff(18:19),'int16')) /32768*10;
                        data(m,12) =  double(typecast(buff(20:21),'int16')) /32768*10;
                        data(m,14) = -double(typecast(buff(22:23),'int16')) /32768*10;
                        data(m,15) =  double(typecast(buff(24:25),'int16')) /10; %temperature
                    case 1
                        data(m,10) = -double(typecast(buff(12:13),'int16')) /50;
                        data(m,9)  = -double(typecast(buff(14:15),'int16')) /50;
                        data(m,11) = -double(typecast(buff(16:17),'int16')) /50;
                        data(m,13) = -double(typecast(buff(18:19),'int16')) /1200;
                        data(m,12) = -double(typecast(buff(20:21),'int16')) /1200;
                        data(m,14) = -double(typecast(buff(22:23),'int16')) /1200;
                        data(m,15) = double(typecast(buff(24:25),'int16')) *0.07386 + 31; %temperature
                end
                data(m,16) = double(buff(26)); %PPS_error
                m = m+1; %指向下一存储位置
                k = k+27; %指向下一帧
            else
                k = k+1; %指向下一字节
            end
        else %到文件尾，退出
            break
        end
    else
        k = k+1;
    end
    if k>n %到文件尾，退出
        break
    end
end

data(m:end,:) = []; %删除空白数据

%% 3.校验cnt，PPS_error
cnt_diff = mod(diff(data(:,1)),256);
if(sum(cnt_diff~=1) ~= 0) %cnt间隔必须是1
    error('cnt error!')
end
if(sum(data(:,16)~=data(1,16)) ~= 0) %PPS_error必须都相同
    error('PPS_error error!')
end

% 统计采样时间
sample_time = mod(diff(data(:,8)),10000); %采样时间，只能为99,100,101，单位：0.1ms
figure
plot(sample_time)
title('采样时间')
sample_time_mean = cumsum(sample_time) ./ (1:length(sample_time))'; %平均采样时间，理论上为10ms，实际可能略高或略低
figure
plot(sample_time_mean)
grid on
title('平均采样时间')

%% 4.提取imu数据，将时间戳转化为GPS时间
n = length(data);
imu_data = zeros(n,7); %IMU数据
% [t, wx,wy,wz, fx,fy,fz], deg/s, g
imu_data(:,2:7) = data(:,9:14);
for k=1:n
    c = data(k,2:7); %时间数组
    c(1) = c(1) + 2000;
    c(4) = c(4) + 8; %转换成北京时间
    if(c(4)>=24) %日期进位
        c(4) = c(4) - 24;
        c(3) = c(3) + 1;
    end
    [gps_week, gps_second] = gps_time(c);
    imu_data(k,1) = gps_second + data(k,8)/10000; %GPS秒
end

% 检查时间是否正确
time_diff = diff(imu_data(:,1));
if(sum(time_diff>0.0102)~=0 || sum(time_diff<0.0098)~=0) %相邻时间应该在10ms
    error('time error!')
end

%% 画imu数据
figure
t = imu_data(:,1) - imu_data(1,1);
subplot(3,2,1)
plot(t,imu_data(:,2))
grid on
set(gca, 'xlim', [t(1),t(end)])
subplot(3,2,3)
plot(t,imu_data(:,3))
grid on
set(gca, 'xlim', [t(1),t(end)])
subplot(3,2,5)
plot(t,imu_data(:,4))
grid on
set(gca, 'xlim', [t(1),t(end)])
subplot(3,2,2)
plot(t,imu_data(:,5))
grid on
set(gca, 'xlim', [t(1),t(end)])
subplot(3,2,4)
plot(t,imu_data(:,6))
grid on
set(gca, 'xlim', [t(1),t(end)])
subplot(3,2,6)
plot(t,imu_data(:,7))
grid on
set(gca, 'xlim', [t(1),t(end)])

%% 5.清除变量
% eval([device,'_data = imu_data;'])
% clearvars imu_data
% clearvars -except *data