clear;
clc;

%% 配置
exe_path = 'F:\GNSS\USRP_B210\Release'; %可执行文件路径
data_path = [cd,'\data']; %数据存储路径，cd返回当前路径
sample_time = 60; %采样时间，s
channel = 1; %通道数量                                                                                       
gain = 38; %增益，dB
ref = 0; %是否使用外部时钟
pps = 0; %是否使用外部PPS，暂时不使用

%% 生成命令
exe_opt = [exe_path,'\gnss_rx'];
data_opt = [' -p "',data_path,'"'];
name_opt = ' -n';
time_opt = [' -t ',num2str(sample_time)];
gain_opt = [' -g ',num2str(gain)];
if channel==1
    channel_opt = [];
else
    channel_opt = ' -d';
end
if ref==0
    ref_opt = [];
else
    ref_opt = ' -r';
end

cmd = [exe_opt, data_opt, name_opt, time_opt, gain_opt, channel_opt, ref_opt]; %控制台命令
system(cmd); %执行系统命令

clearvars