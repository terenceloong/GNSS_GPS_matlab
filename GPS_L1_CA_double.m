% GPS双天线定位测姿程序
% 时钟、频率反馈，在准确GPS时间点进行定位、测姿

clear
clc

%% 计时开始
tic

%% 创建日志文件
fclose('all'); %关闭之前打开的所有文件
logID_A = fopen('logA.txt', 'w'); %创建日志文件（时间顺序的日志）
logID_B = fopen('logB.txt', 'w');

%% 文件路径
file_path_A = '.\data\7_2\data_20190702_111609_ch1.dat';
file_path_B = [file_path_A(1:(end-5)),'2.dat'];
plot_gnss_file(file_path_A); %显示前0.1s数据
plot_gnss_file(file_path_B);
sample_offset = 0*4e6; %抛弃前多少个采样点

%% 运行时间
msToProcess = 60*1000; %处理总时间
sampleFreq = 4e6; %接收机采样频率

%% 参考位置
p0 = [45.730952, 126.624970, 212]; %2A楼顶

%% 数据缓存
buffBlkNum = 40;                     %采样数据缓存块数量（要保证捕获时存储恰好从头开始）
buffBlkSize = 4000;                  %一个块的采样点数（1ms）
buffSize = buffBlkSize * buffBlkNum; %采样数据缓存大小
buff_A = zeros(2,buffSize);          %采样数据缓存，第一行I，第二行Q
buff_B = zeros(2,buffSize);
buffBlkPoint = 0;                    %数据该往第几块存，从0开始
buffHead = 0;                        %最新数据的序号，buffBlkSize的倍数

%% 获取文件时间
tf = sscanf(file_path_A((end-22):(end-8)), '%4d%02d%02d_%02d%02d%02d')'; %数据文件开始采样时间（日期时间数组）
[tw, ts] = gps_time(tf); %tw：GPS周数，ts：GPS周内秒数
ta = [ts,0,0] + sample2dt(sample_offset, sampleFreq); %初始化接收机时间，[s,ms,us]
ta = time_carry(round(ta,2)); %取整

%% 根据历书获取当前可能见到的卫星（*）
% svList = [2;6;12];
svList = gps_constellation(tf, p0);
svN = length(svList);

%% 为每颗可能见到的卫星分配跟踪通道
channels_A = repmat(GPS_L1_CA_channel_struct(), svN,1);
channels_B = repmat(GPS_L1_CA_channel_struct(), svN,1);
for k=1:svN
    channels_A(k).PRN = svList(k);
    channels_A(k).state = 0; %状态未激活
    channels_B(k).PRN = svList(k);
    channels_B(k).state = 0; %状态未激活
end

%% 预置星历
ephemeris_file = ['./ephemeris/',file_path_A((end-22):(end-8)),'.mat'];
if exist(ephemeris_file, 'file')
    load(ephemeris_file); %星历存在，加载星历文件，星历变量名为ephemeris，星历为列；电离层参数变量名为ion
else
    ephemeris = NaN(26,32); %星历不存在，设置空的星历
    ion = NaN(1,8); %空的电离层参数
end
for k=1:svN
    PRN = svList(k);
    channels_A(k).ephemeris = ephemeris(:,PRN); %为通道的星历赋值
    channels_B(k).ephemeris = ephemeris(:,PRN);
    if ~isnan(ephemeris(1,PRN)) %如果存在某颗卫星的星历，打印日志
        fprintf(logID_A, '%2d: Load ephemeris.\r\n', PRN);
        fprintf(logID_B, '%2d: Load ephemeris.\r\n', PRN);
    end
end

%% 创建跟踪结果存储空间
% 分配了msToProcess行，每跟踪一次输出一次结果，最后删除多余的行
trackResults_A = repmat(trackResult_struct(msToProcess), svN,1);
trackResults_B = repmat(trackResult_struct(msToProcess), svN,1);
for k=1:svN
    trackResults_A(k).PRN = svList(k);
    trackResults_B(k).PRN = svList(k);
end

%% 接收机状态
receiverState = 0; %接收机状态，0表示未初始化，时间还不对，1表示时间已经校正
deltaFreq = 0; %时钟差，理解为百分比，如果差1e-9，生成1500e6Hz的波会差1.5Hz
dtpos = 10; %定位时间间隔，ms
tp = [ta(1),0,0]; %tp为下次定位时间
tp(2) = (floor(ta(2)/dtpos)+1) * dtpos; %加到下个整目标时间
tp = time_carry(tp); %进位

%% 创建接收机输出存储空间
% 分配msToProcess/dtpos行，接收机初始化后每dtpos ms输出一次，最后删除多余的行
output_ta = zeros(msToProcess/dtpos,1); %时间，ms
output_pos = zeros(msToProcess/dtpos,8); %定位，[位置、速度、钟差、钟频差]
output_sv = zeros(svN,8,msToProcess/dtpos); %卫星信息，[位置、伪距、速度、伪距率]
output_df = zeros(msToProcess/dtpos,1); %修正用的钟频差（滤波后的钟频差）
output_dphase = NaN(msToProcess/dtpos,svN); %相位差
no = 1; %指向当前存储行

%% 打开文件，创建进度条
% 文件A
fileID_A = fopen(file_path_A, 'r');
fseek(fileID_A, round(sample_offset*4), 'bof');
if int64(ftell(fileID_A))~=int64(sample_offset*4)
    error('Sample offset error!');
end
% 文件B
fileID_B = fopen(file_path_B, 'r');
fseek(fileID_B, round(sample_offset*4), 'bof');
if int64(ftell(fileID_B))~=int64(sample_offset*4)
    error('Sample offset error!');
end
% 进度条
f = waitbar(0, ['0s/',num2str(msToProcess/1000),'s']);

%% 信号处理
for t=1:msToProcess
    % 更新进度条
    if mod(t,1000)==0 %1s步进
        waitbar(t/msToProcess, f, [num2str(t/1000),'s/',num2str(msToProcess/1000),'s']);
    end
    
    % 读数据
    buff_A(:,buffBlkPoint*buffBlkSize+(1:buffBlkSize)) = double(fread(fileID_A, [2,buffBlkSize], 'int16')); %天线A
    buff_B(:,buffBlkPoint*buffBlkSize+(1:buffBlkSize)) = double(fread(fileID_B, [2,buffBlkSize], 'int16')); %天线B
    buffBlkPoint = buffBlkPoint + 1;
    buffHead = buffBlkPoint * buffBlkSize;
    if buffBlkPoint==buffBlkNum
        buffBlkPoint = 0; %缓存从头开始
    end
    
	%% 更新接收机时间
    % 当前最后一个采样的接收机时间
    sampleFreq_real = sampleFreq * (1+deltaFreq); %真实的采样频率
    ta = time_carry(ta + sample2dt(buffBlkSize, sampleFreq_real));
    
    %% 捕获
    % 每1s的采样点搜索一次
    if mod(t,1000)==0
        for k=1:svN %搜索所有可能见到的卫星
            %====天线A
            if channels_A(k).state==0 %如果通道未激活，捕获尝试激活
                [acqResult, peakRatio] = GPS_L1_CA_acq_one(svList(k), buff_A(:,(end-2*8000+1):end)); %2ms数据捕获
                if ~isempty(acqResult) %成功捕获
                    channels_A(k) = GPS_L1_CA_channel_init(channels_A(k), acqResult, t*buffBlkSize, sampleFreq); %激活通道
                    fprintf(logID_A, '%2d: Acquired at %ds, peakRatio=%.2f\r\n', svList(k), t/1000, peakRatio); %打印捕获日志
                end
            end
            %====天线B
            if channels_B(k).state==0 %如果通道未激活，捕获尝试激活
                [acqResult, peakRatio] = GPS_L1_CA_acq_one(svList(k), buff_B(:,(end-2*8000+1):end)); %2ms数据捕获
                if ~isempty(acqResult) %成功捕获
                    channels_B(k) = GPS_L1_CA_channel_init(channels_B(k), acqResult, t*buffBlkSize, sampleFreq); %激活通道
                    fprintf(logID_B, '%2d: Acquired at %ds, peakRatio=%.2f\r\n', svList(k), t/1000, peakRatio); %打印捕获日志
                end
            end
        end
    end
    
    %% 跟踪
    for k=1:svN
        %====天线A
        if channels_A(k).state~=0 %如果通道激活，进行跟踪
            while 1
                % 判断是否有完整的跟踪数据
                if mod(buffHead-channels_A(k).trackDataHead,buffSize)>(buffSize/2)
                    break
                end
                % 存跟踪结果（通道参数）
                n = trackResults_A(k).n;
                trackResults_A(k).dataIndex(n,:)    = channels_A(k).dataIndex;
                trackResults_A(k).ts0(n,:)          = channels_A(k).ts0;
                trackResults_A(k).remCodePhase(n,:) = channels_A(k).remCodePhase;
                trackResults_A(k).codeFreq(n,:)     = channels_A(k).codeFreq;
                trackResults_A(k).remCarrPhase(n,:) = channels_A(k).remCarrPhase;
                trackResults_A(k).carrFreq(n,:)     = channels_A(k).carrFreq;
                % 基带处理
                trackDataHead = channels_A(k).trackDataHead;
                trackDataTail = channels_A(k).trackDataTail;
                if trackDataHead>trackDataTail
                    [channels_A(k), I_Q, disc, bitStartFlag] = ...
                        GPS_L1_CA_track(channels_A(k), sampleFreq_real, buffSize, buff_A(:,trackDataTail:trackDataHead), logID_A);
                else
                    [channels_A(k), I_Q, disc, bitStartFlag] = ...
                        GPS_L1_CA_track(channels_A(k), sampleFreq_real, buffSize, [buff_A(:,trackDataTail:end),buff_A(:,1:trackDataHead)], logID_A);
                end
                % 存跟踪结果（跟踪结果）
                trackResults_A(k).I_Q(n,:)          = I_Q;
                trackResults_A(k).disc(n,:)         = disc;
                trackResults_A(k).bitStartFlag(n,:) = bitStartFlag;
                trackResults_A(k).CN0(n,:)          = channels_A(k).CN0;
                trackResults_A(k).carrAcc(n,:)      = channels_A(k).carrAcc;
                trackResults_A(k).Px(n,:)           = sqrt(diag(channels_A(k).Px)')*3;
                trackResults_A(k).n                 = n + 1;
            end
        end
        %====天线B
        if channels_B(k).state~=0 %如果通道激活，进行跟踪
            while 1
                % 判断是否有完整的跟踪数据
                if mod(buffHead-channels_B(k).trackDataHead,buffSize)>(buffSize/2)
                    break
                end
                % 存跟踪结果（通道参数）
                n = trackResults_B(k).n;
                trackResults_B(k).dataIndex(n,:)    = channels_B(k).dataIndex;
                trackResults_B(k).ts0(n,:)          = channels_B(k).ts0;
                trackResults_B(k).remCodePhase(n,:) = channels_B(k).remCodePhase;
                trackResults_B(k).codeFreq(n,:)     = channels_B(k).codeFreq;
                trackResults_B(k).remCarrPhase(n,:) = channels_B(k).remCarrPhase;
                trackResults_B(k).carrFreq(n,:)     = channels_B(k).carrFreq;
                % 基带处理
                trackDataHead = channels_B(k).trackDataHead;
                trackDataTail = channels_B(k).trackDataTail;
                if trackDataHead>trackDataTail
                    [channels_B(k), I_Q, disc, bitStartFlag] = ...
                        GPS_L1_CA_track(channels_B(k), sampleFreq_real, buffSize, buff_B(:,trackDataTail:trackDataHead), logID_B);
                else
                    [channels_B(k), I_Q, disc, bitStartFlag] = ...
                        GPS_L1_CA_track(channels_B(k), sampleFreq_real, buffSize, [buff_B(:,trackDataTail:end),buff_B(:,1:trackDataHead)], logID_B);
                end
                % 存跟踪结果（跟踪结果）
                trackResults_B(k).I_Q(n,:)          = I_Q;
                trackResults_B(k).disc(n,:)         = disc;
                trackResults_B(k).bitStartFlag(n,:) = bitStartFlag;
                trackResults_B(k).CN0(n,:)          = channels_B(k).CN0;
                trackResults_B(k).carrAcc(n,:)      = channels_B(k).carrAcc;
                trackResults_B(k).Px(n,:)           = sqrt(diag(channels_B(k).Px)')*3;
                trackResults_B(k).n                 = n + 1;
            end
        end
    end
    
    %% 检查目标时间是否到达
    dtp = (ta(1)-tp(1)) + (ta(2)-tp(2))/1e3 + (ta(3)-tp(3))/1e6; %当前采样时间与定位时间之差，>=0时表示当前采样时间已经到达或超过定位时间
    
    %% 定位
    % 只使用A天线
    if dtp>=0
        %--------计算卫星信息
        sv = NaN(svN,8);
        for k=1:svN
            if channels_A(k).state==2 %检查所有通道状态，对跟踪到的通道计算卫星信息，[位置、伪距、速度、伪距率]
                dn = mod(buffHead-channels_A(k).trackDataTail+1, buffSize) - 1; %trackDataTail恰好超前buffHead一个时，dn=-1
                dtc = dn / sampleFreq_real; %当前采样时间与跟踪点的时间差
                carrFreq = channels_A(k).carrFreq + 1575.42e6*deltaFreq; %修正后的载波频率
                codeFreq = (carrFreq/1575.42e6+1)*1.023e6; %通过载波频率计算的码频率
                codePhase = channels_A(k).remCodePhase + (dtc-dtp)*codeFreq; %定位点码相位
                ts0 = [floor(channels_A(k).ts0/1e3), mod(channels_A(k).ts0,1e3), 0] + [0, floor(codePhase/1023), mod(codePhase/1023,1)*1e3]; %定位点的码发射时间
                [sv(k,:),~] = sv_ecef(channels_A(k).ephemeris, tp, ts0); %根据星历计算卫星[位置、伪距、速度]
                sv(k,8) = -carrFreq/1575.42e6*299792458;%载波频率转化为速度
            end
        end
        %--------定位
        sv_visible = sv(~isnan(sv(:,1)),:); %提取可见卫星
        pos = pos_solve(sv_visible); %定位，如果不够4颗卫星返回8个NaN
        %--------检查初始化
        if receiverState==0
            if ~isnan(pos(7)) %钟差不为NaN
                if abs(pos(7))>0.1e-3 %钟差大于0.1ms，修正接收机时间
                    ta = ta - sec2smu(pos(7)); %时钟修正
                    ta = time_carry(ta);
                    tp(1) = ta(1); %更新下次定位时间
                    tp(2) = (floor(ta(2)/dtpos)+1) * dtpos;
                    tp = time_carry(tp);
                else %钟差小于0.1ms，初始化结束
                    receiverState = 1;
                end
            end
        end
        %--------接收机已完成初始化
        if receiverState==1
            %--------时钟反馈修正
            if ~isnan(pos(7)) %钟差不为NaN
                deltaFreq = deltaFreq + 10*pos(8)*dtpos/1000; %钟频差累加
                ta = ta - 10*sec2smu(pos(7))*dtpos/1000; %时钟修正（可以不用进位，在下次更新时进位）
            end
            %--------存储输出
            output_ta(no) = tp(1)*1000 + tp(2); %时间戳，ms
            output_pos(no,:) = pos;
            output_sv(:,:,no) = sv;
            output_df(no) = deltaFreq;
            % no = no + 1;
        end
    end
    
    %% 计算相位差
    % A相位 - B相位
    if dtp>=0 && receiverState==1
        for k=1:svN
            if channels_A(k).state==2 && channels_B(k).state==2 %两个天线都跟踪到该颗卫星
                % 天线A
                dn = mod(buffHead-channels_A(k).trackDataTail+1, buffSize) - 1;
                dtc = dn / sampleFreq_real;
                dt = dtc - dtp;
                phase_A = channels_A(k).remCarrPhase + channels_A(k).carrFreq*dt + 0.5*channels_A(k).carrAcc*dt^2; %载波相位
                % 天线B
                dn = mod(buffHead-channels_B(k).trackDataTail+1, buffSize) - 1;
                dtc = dn / sampleFreq_real;
                dt = dtc - dtp;
                phase_B = channels_B(k).remCarrPhase + channels_B(k).carrFreq*dt + 0.5*channels_B(k).carrAcc*dt^2; %载波相位
                % 相位差
                if channels_A(k).inverseFlag*channels_B(k).inverseFlag==1 %两个天线相位翻转相同
                    output_dphase(no,k) = mod((channels_A(k).carrCirc+phase_A)-(channels_B(k).carrCirc+phase_B)    +500,1000) - 500;
                else %两个天线相位翻转不同
                    output_dphase(no,k) = mod((channels_A(k).carrCirc+phase_A)-(channels_B(k).carrCirc+phase_B)+0.5+500,1000) - 500;
                end
            end
        end
        no = no + 1;
    end
    
    %% 更新下次目标时间
    if dtp>=0
        tp = time_carry(tp + [0,dtpos,0]);
    end
    
end

%% 关闭文件，关闭进度条
fclose(fileID_A);
fclose(fileID_B);
fclose(logID_A);
fclose(logID_B);
close(f);

%% 删除空白数据
for k=1:svN
    trackResults_A(k) = trackResult_clean(trackResults_A(k));
    trackResults_B(k) = trackResult_clean(trackResults_B(k));
end
output_ta(no:end) = [];
output_pos(no:end,:) = [];
output_sv(:,:,no:end) = [];
output_df(no:end) = [];
output_dphase(no:end,:) = [];

%% 打印通道日志（*）
clc
disp('<--------antenna A-------->')
print_log('logA.txt', svList);
disp('<--------antenna B-------->')
print_log('logB.txt', svList);

%% 保存星历
% 每次运行完都会保存，有新星历自动添加
for k=1:svN
    PRN = channels_A(k).PRN;
    if isnan(ephemeris(1,PRN)) %星历文件中没有
        if ~isnan(channels_A(k).ephemeris(1))
            ephemeris(:,PRN) = channels_A(k).ephemeris; %插入星历
        elseif ~isnan(channels_B(k).ephemeris(1))
            ephemeris(:,PRN) = channels_B(k).ephemeris; %插入星历
        end
    end
end
save(ephemeris_file, 'ephemeris', 'ion');

%% 画图（*）
for k=1:svN
    if trackResults_A(k).n==1 && trackResults_B(k).n==1 %不画没跟踪的通道
        continue
    end
    
    % 建立坐标轴
    screenSize = get(0,'ScreenSize'); %获取屏幕尺寸
    if screenSize(3)==1920 %根据屏幕尺寸设置画图范围
        figure('Position', [390, 280, 1140, 670]);
    elseif screenSize(3)==1368
        figure('Position', [114, 100, 1140, 670]);
    else
        error('Screen size error!')
    end
    ax1 = axes('Position', [0.08, 0.4, 0.38, 0.53]);
    hold(ax1,'on');
    axis(ax1, 'equal');
    title(['PRN = ',num2str(svList(k))])
    ax2 = axes('Position', [0.53, 0.7 , 0.42, 0.25]);
    hold(ax2,'on');
    ax3 = axes('Position', [0.53, 0.38, 0.42, 0.25]);
    hold(ax3,'on');
    ax4 = axes('Position', [0.53, 0.06, 0.42, 0.25]);
    hold(ax4,'on');
    grid(ax4,'on');
    ax5 = axes('Position', [0.05, 0.06, 0.42, 0.25]);
    hold(ax5,'on');
    grid(ax5,'on');
    
    % 画图
    plot(ax1, trackResults_A(k).I_Q(1001:end,1),trackResults_A(k).I_Q(1001:end,4), 'LineStyle','none', 'Marker','.', 'Color',[0,0.447,0.741])
    plot(ax2, trackResults_A(k).dataIndex/sampleFreq, trackResults_A(k).I_Q(:,1), 'Color',[0,0.447,0.741])
    
%     index = find(trackResults_A(k).bitStartFlag==double('H')); %寻找帧头阶段（粉色）
%     plot(ax2, trackResults_A(k).dataIndex(index)/sampleFreq, trackResults_A(k).I_Q(index,1), 'LineStyle','none', 'Marker','.', 'Color','m')
%     index = find(trackResults_A(k).bitStartFlag==double('C')); %校验帧头阶段（蓝色）
%     plot(ax2, trackResults_A(k).dataIndex(index)/sampleFreq, trackResults_A(k).I_Q(index,1), 'LineStyle','none', 'Marker','.', 'Color','b')
%     index = find(trackResults_A(k).bitStartFlag==double('E')); %解析星历阶段（红色）
%     plot(ax2, trackResults_A(k).dataIndex(index)/sampleFreq, trackResults_A(k).I_Q(index,1), 'LineStyle','none', 'Marker','.', 'Color','r')
    %---------------------------------------------------------------------%
    plot(ax1, trackResults_B(k).I_Q(1001:end,1),trackResults_B(k).I_Q(1001:end,4), 'LineStyle','none', 'Marker','.', 'Color',[0.850,0.325,0.098])
    plot(ax3, trackResults_B(k).dataIndex/sampleFreq, trackResults_B(k).I_Q(:,1), 'Color',[0.850,0.325,0.098])
    
%     index = find(trackResults_B(k).bitStartFlag==double('H')); %寻找帧头阶段（粉色）
%     plot(ax3, trackResults_B(k).dataIndex(index)/sampleFreq, trackResults_B(k).I_Q(index,1), 'LineStyle','none', 'Marker','.', 'Color','m')
%     index = find(trackResults_B(k).bitStartFlag==double('C')); %校验帧头阶段（蓝色）
%     plot(ax3, trackResults_B(k).dataIndex(index)/sampleFreq, trackResults_B(k).I_Q(index,1), 'LineStyle','none', 'Marker','.', 'Color','b')
%     index = find(trackResults_B(k).bitStartFlag==double('E')); %解析星历阶段（红色）
%     plot(ax3, trackResults_B(k).dataIndex(index)/sampleFreq, trackResults_B(k).I_Q(index,1), 'LineStyle','none', 'Marker','.', 'Color','r')

    plot(ax4, trackResults_A(k).dataIndex/sampleFreq, trackResults_A(k).carrFreq, 'LineWidth',1.5, 'Color',[0,0.447,0.741]) %载波频率
    plot(ax4, trackResults_B(k).dataIndex/sampleFreq, trackResults_B(k).carrFreq, 'LineWidth',1.5, 'Color',[0.850,0.325,0.098])
    
    plot(ax5, trackResults_A(k).dataIndex/sampleFreq, trackResults_A(k).carrAcc, 'Color',[0,0.447,0.741]) %视线方向加速度
    plot(ax5, trackResults_B(k).dataIndex/sampleFreq, trackResults_B(k).carrAcc, 'Color',[0.850,0.325,0.098])
    
    % 调整坐标轴
    set(ax2, 'XLim',[0,msToProcess/1000])
    set(ax3, 'XLim',[0,msToProcess/1000])

    ax2_ylim = get(ax2, 'YLim');
    ax3_ylim = get(ax3, 'YLim');
    ylim = max(abs([ax2_ylim,ax3_ylim]));
    set(ax2, 'YLim',[-ylim,ylim])
    set(ax3, 'YLim',[-ylim,ylim])
    
    set(ax4, 'XLim',[0,msToProcess/1000])
    set(ax5, 'XLim',[0,msToProcess/1000])
end

clearvars k screenSize ax1 ax2 ax3 ax4 ax5 index ax2_ylim ax3_ylim ylim

%% 画相位差（*）
colorTable = [    0, 0.447, 0.741;
              0.850, 0.325, 0.098;
              0.929, 0.694, 0.125;
              0.494, 0.184, 0.556;
              0.466, 0.674, 0.188;
              0.301, 0.745, 0.933;
              0.635, 0.078, 0.184;
                  0,     0,     1;
                  1,     0,     0;
                  0,     1,     0;
                  0,     0,     0;];
figure
hold on
grid on
legend_str = [];
for k=1:svN
    if sum(~isnan(output_dphase(:,k)))~=0
        plot(output_dphase(:,k), 'LineWidth',1, 'Color',colorTable(k,:))
        eval('legend_str = [legend_str; string(num2str(svList(k)))];')
    end
end
legend(legend_str)
% set(gca,'Ylim',[0,1])
title('Phase difference')
clearvars colorTable legend_str k

%% 清除变量（*）
clearvars -except sampleFreq msToProcess ...
                  p0 tf svList svN ...
                  channels_A trackResults_A ...
                  channels_B trackResults_B ...
                  output_ta output_pos output_sv output_df output_dphase ...
                  ion
              
save result_double.mat

%% 计时结束
toc