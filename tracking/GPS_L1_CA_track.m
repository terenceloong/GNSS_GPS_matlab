function [channel, I_Q, disc, bitStartFlag] = GPS_L1_CA_track(channel, sampleFreq, buffSize, rawSignal, logID)

bitStartFlag = 0;

%% 提取通道信息（跟踪算法用到的控制参数）
trackStage     = channel.trackStage;
msgStage       = channel.msgStage;
cnt_t          = channel.cnt_t;
cnt_m          = channel.cnt_m;
code           = channel.code;
timeIntMs      = channel.timeIntMs;
blkSize        = channel.blkSize;
carrNco        = channel.carrNco;
codeNco        = channel.codeNco;
carrAcc        = channel.carrAcc;
remCarrPhase   = channel.remCarrPhase;
remCodePhase   = channel.remCodePhase;
carrCirc       = channel.carrCirc;
I_P0           = channel.I_P0;
Q_P0           = channel.Q_P0;
FLL            = channel.FLL;
PLL            = channel.PLL;
DLL            = channel.DLL;
bitSyncTable   = channel.bitSyncTable;
bitBuff        = channel.bitBuff;
frameBuff      = channel.frameBuff;
frameBuffPoint = channel.frameBuffPoint;
CN0            = channel.CN0;
P              = channel.Px;

channel.dataIndex = channel.dataIndex + blkSize; %下个数据段开始的数据索引
channel.ts0       = channel.ts0 + timeIntMs; %当前码周期结束的时间，ms

timeInt = timeIntMs * 0.001; %积分时间，s
codeInt = timeIntMs * 1023; %积分时间内码片个数
pointInt = 20 / timeIntMs; %一个电文比特中有多少个积分点

%% 基本处理
% 时间序列
t = (0:blkSize-1) / sampleFreq;
te = blkSize / sampleFreq;

% 生成本地载波
theta = (remCarrPhase+(carrNco+0.5*carrAcc*t).*t) * 2; %变频率载波，乘2因为后面是以pi为单位求三角函数
carr_cos = cospi(theta); %本地载波
carr_sin = sinpi(theta);
theta_next = remCarrPhase + carrNco*te + 0.5*carrAcc*te^2;
remCarrPhase = mod(theta_next, 1); %剩余载波相位，周
carrCirc = mod(floor(carrCirc+theta_next), 1000); %载波经过的整周数

% 生成本地码
tcode = remCodePhase + codeNco*t + 2; %加2保证求滞后码时大于1
earlyCode  = code(floor(tcode+0.5)); %超前码
promptCode = code(floor(tcode));     %即时码
lateCode   = code(floor(tcode-0.5)); %滞后码
remCodePhase = remCodePhase + codeNco*te - codeInt; %剩余载波相位，周

% 原始数据乘载波
% 将复数相乘拆开，为了避免取实部、虚部时调用函数耗时
iBasebandSignal = rawSignal(1,:).*carr_cos + rawSignal(2,:).*carr_sin; %乘负载波
qBasebandSignal = rawSignal(2,:).*carr_cos - rawSignal(1,:).*carr_sin;

% 六路积分
% 用行向量与列向量相乘代替sum(X.*X)，加速
I_E = iBasebandSignal * earlyCode;
Q_E = qBasebandSignal * earlyCode;
I_P = iBasebandSignal * promptCode;
Q_P = qBasebandSignal * promptCode;
I_L = iBasebandSignal * lateCode;
Q_L = qBasebandSignal * lateCode;

% 码鉴相器
S_E = sqrt(I_E^2+Q_E^2);
S_L = sqrt(I_L^2+Q_L^2);
codeError = 0.5 * (S_E-S_L)/(S_E+S_L); %单位：码片
[channel.codeStd, codeSigma] = std_rec(channel.codeStd ,codeError); %计算码鉴相器误差标准差

% 载波鉴相器
carrError = atan(Q_P/I_P) / (2*pi); %单位：周
[channel.carrStd, carrSigma] = std_rec(channel.carrStd ,carrError); %计算载波鉴相器误差标准差

% 鉴频器
if ~isnan(I_P0)
    yc = I_P0*I_P + Q_P0*Q_P;
    ys = I_P0*Q_P - Q_P0*I_P;
    freqError = atan(ys/yc)/timeInt / (2*pi); %单位：Hz
else
    freqError = 0;
end

%% 跟踪算法
switch trackStage
    case 'F' %<<====频率牵引
        %----FLL
        FLL.Int = FLL.Int + FLL.K*freqError*timeInt; %锁频环积分器
        carrNco = FLL.Int;
        carrFreq = FLL.Int;
        % 500ms后转到传统跟踪
        cnt_t = cnt_t + 1;
        if cnt_t==500
            cnt_t = 0; %计数器清零
            PLL.Int = FLL.Int; %初始化锁相环积分器
            trackStage = 'T';
            fprintf(logID, '%2d: Start traditional tracking at %.8fs\r\n', ...
                    channel.PRN, channel.dataIndex/sampleFreq);
        end
        %----DLL
        DLL.Int = DLL.Int + DLL.K2*codeError*timeInt; %延迟锁定环积分器
        codeNco = DLL.Int + DLL.K1*codeError;
        codeFreq = DLL.Int;
        
	case 'T' %<<====传统跟踪
        %----PLL
        PLL.Int = PLL.Int + PLL.K2*carrError*timeInt; %锁相环积分器
        carrNco = PLL.Int + PLL.K1*carrError;
        carrFreq = PLL.Int;
        % 500ms后进行比特同步
        if msgStage=='I'
            cnt_t = cnt_t + 1;
            if cnt_t==500
                cnt_t = 0; %计数器清零
                msgStage = 'B';
                fprintf(logID, '%2d: Start bit synchronization at %.8fs\r\n', ...
                        channel.PRN, channel.dataIndex/sampleFreq);
            end
        end
        %----DLL
        DLL.Int = DLL.Int + DLL.K2*codeError*timeInt; %延迟锁定环积分器
        codeNco = DLL.Int + DLL.K1*codeError;
        codeFreq = DLL.Int;
        
    case 'K' %<<====卡尔曼滤波跟踪
        Phi = eye(4);
        Phi(1,3) = timeInt/1540;
        Phi(1,4) = 0.5*timeInt^2/1540;
        Phi(2,3) = timeInt;
        Phi(2,4) = 0.5*timeInt^2;
        Phi(3,4) = timeInt;
        H = zeros(2,4);
        H(1,1) = 1;
        H(2,2) = 1;
        H(2,3) = -timeInt/2;
        % Q阵的两个参数需要调节，第一个对应钟频的噪声，第二个对应钟漂移的噪声
%         Q = diag([0, 0, 3, 0.02].^2) * timeInt^2;
%         Q = diag([0, 0, 1, 0.01].^2) * timeInt^2;
        Q = diag([0, 0, 4, 2].^2) * timeInt^2;
        R = diag([codeSigma, carrSigma].^2);
        Z = [codeError; carrError];
        % 卡尔曼滤波
        P = Phi*P*Phi' + Q;
        K = P*H' / (H*P*H'+R);
        X = K*Z;
        P = (eye(4)-K*H)*P;
        P = (P+P')/2;
        % 修正
        remCodePhase = remCodePhase + X(1);
        remCarrPhase = remCarrPhase + X(2);
        carrNco = carrNco + carrAcc*(blkSize/sampleFreq) + X(3);
        codeNco = 1.023e6 + carrNco/1540;
        carrAcc = carrAcc + X(4);
        carrFreq = carrNco;
        codeFreq = codeNco;
        
	otherwise
end

%% 电文解析算法
switch msgStage %I, B, W, H, C, E
    case 'I' %<<====空闲
        
    case 'B' %<<====比特同步
        % 必须是1ms积分时间，进行2s，有100个比特
        % 比特同步后可以实现更长的积分时间
        % 比特同步后可以进行载噪比计算
        cnt_m = cnt_m + 1;
        if (I_P0*I_P)<0 %发现电平翻转
            index = mod(cnt_m-1,20) + 1;
            bitSyncTable(index) = bitSyncTable(index) + 1; %统计表中的对应位加1
        end
        if cnt_m==2000 %2s后检验统计表
            if max(bitSyncTable)>15 && (sum(bitSyncTable)-max(bitSyncTable))<=5 %确定电平翻转位置（电平翻转大都发生在一个点上）
                %------------------------------------------------------------------------------------%
%                 trackStage = 'K'; %比特同步后转到卡尔曼滤波跟踪
%                 fprintf(logID, '%2d: Start Kalman tracking at %.8fs\r\n', ...
%                         channel.PRN, channel.dataIndex/sampleFreq);
                %------------------------------------------------------------------------------------%
                [~,cnt_m] = max(bitSyncTable);
                bitSyncTable = zeros(1,20); %比特同步统计表清零
                cnt_m = -cnt_m + 1;
                if cnt_m==0
                    msgStage = 'H'; %进入寻找帧头模式
                    fprintf(logID, '%2d: Start find head at %.8fs\r\n', ...
                            channel.PRN, channel.dataIndex/sampleFreq);
                else
                    msgStage = 'W'; %进入等待cnt_m==0模式
                end
            else
                channel.state = 0; %比特同步失败，关闭通道
                fprintf(logID, '%2d: ***Bit synchronization fails at %.8fs\r\n', ...
                        channel.PRN, channel.dataIndex/sampleFreq);
            end
        end
        
    case 'W' %<<====等待cnt_m==0
        cnt_m = cnt_m + 1;
        if cnt_m==0
            msgStage = 'H'; %进入寻找帧头模式
            fprintf(logID, '%2d: Start find head at %.8fs\r\n', ...
                    channel.PRN, channel.dataIndex/sampleFreq);
        end
        
    otherwise %<<====已经完成比特同步
        cnt_m = cnt_m + 1;
        bitBuff(1,cnt_m) = I_P; %往比特缓存中存数
        bitBuff(2,cnt_m) = Q_P; %往比特缓存中存数
        if cnt_m==1 %标记当前跟踪的数据段为比特开始位置
            bitStartFlag = double(msgStage);
        end
        if cnt_m==pointInt %跟踪完一个比特
            cnt_m = 0; %计数器清零
            %-------------------------------------------------------------%
            % 计算载噪比
            Ps = bitBuff(1,1:pointInt).^2 + bitBuff(2,1:pointInt).^2; %每个点的功率
            WBP = sum(Ps); %宽带功率，所有点的功率求和（先平方再求和）
            Is = sum(bitBuff(1,1:pointInt)); %合成I
            Qs = sum(bitBuff(2,1:pointInt)); %合成Q
            NBP = Is^2 + Qs^2; %窄带功率，合成IQ的功率，信号越好，窄带功率越大（先求和再平方）
            if CN0==0 %初始化求均值结构体
                channel.NWmean.buff = ones(1,channel.NWmean.buffSize)*(NBP/WBP);
                channel.NWmean.E0 = NBP/WBP;
            end
            [channel.NWmean, NWm] = mean_rec(channel.NWmean ,NBP/WBP); %计算Z的均值
            CN0 = 10*log10((NWm-1)/(pointInt-NWm)/timeInt); %载噪比
            if CN0<=35 %判定为失锁
                channel.state = 0;
                fprintf(logID, '%2d: ***Lose lock at %.8fs\r\n', ...
                        channel.PRN, channel.dataIndex/sampleFreq);
            end
            %-------------------------------------------------------------%
            bit = sum(bitBuff(1,1:pointInt)) > 0; %判断比特值，0/1
            frameBuffPoint = frameBuffPoint + 1;
            frameBuff(frameBuffPoint) = (double(bit) - 0.5) * 2; %存储比特值，±1
            switch msgStage
                case 'H' %<<====寻找帧头
                    if frameBuffPoint>=10 %至少有10个比特，前两个用来校验
                        if abs(sum(frameBuff(frameBuffPoint+(-7:0)).*[1,-1,-1,-1,1,-1,1,1]))==8 %检测到疑似帧头
                            frameBuff(1:10) = frameBuff(frameBuffPoint+(-9:0)); %将帧头提前
                            frameBuffPoint = 10;
                            msgStage = 'C'; %进入校验帧头模式
                        end
                    end
                    if frameBuffPoint==1502 %防止Bug，一般到不了这里，30s还没找到帧头早就被判定为失锁了
                        frameBuffPoint = 0;
                    end
                case 'C' %<<====校验帧头
                    if frameBuffPoint==310 %存储了一个子帧，2+300+8
                        if GPS_L1_CA_check(frameBuff(1:32))==1 && GPS_L1_CA_check(frameBuff(31:62))==1  && ... %校验通过
                            abs(sum(frameBuff(303:310).*[1,-1,-1,-1,1,-1,1,1]))==8
                            % 获取电文时间
                            % frameBuff(32)为上一字的最后一位，校验时控制电平翻转，为1表示翻转，为0表示不翻转，参见ICD-GPS最后几页
                            bits = -frameBuff(32) * frameBuff(33:49); %电平翻转，31~47比特
                            bits = dec2bin(bits>0)'; %±1数组转化为01字符串
                            TOW = bin2dec(bits); %01字符串转换为十进制数
                            channel.ts0 = (TOW*6+0.16)*1000; %ms，0.16=8/50
                            channel.inverseFlag = frameBuff(62); %相位翻转标志，1表示翻转，-1表示不翻转
                            if ~isnan(channel.ephemeris(1))
                                channel.state = 2; %更新状态（知道码发射时间，而且有星历）
                            end
                            msgStage = 'E'; %进入解析星历模式
                            fprintf(logID, '%2d: Start parse ephemeris at %.8fs\r\n', ...
                                    channel.PRN, channel.dataIndex/sampleFreq);
                        else %校验未通过
                            for k=11:310 %检查其他比特中有没有帧头
                                if abs(sum(frameBuff(k+(-7:0)).*[1,-1,-1,-1,1,-1,1,1]))==8 %检测到疑似帧头
                                    frameBuff(1:320-k) = frameBuff(k-9:310); %将帧头和后面的比特提前，320-k = 310-(k-9)+1
                                    frameBuffPoint = 320-k;
                                    break
                                end
                            end
                            if frameBuffPoint==310 %没检测到疑似帧头
                                frameBuff(1:9) = frameBuff(302:310); %将未检测的比特提前
                                frameBuffPoint = 9;
                                msgStage = 'H'; %再次寻找帧头
                            end
                        end
                    end
                case 'E' %<<====解析星历
                    if frameBuffPoint==1502 %跟踪完5帧
                        ephemeris = GPS_L1_CA_ephemeris(frameBuff); %解析星历
                        if isempty(ephemeris) %星历错误
                            fprintf(logID, '%2d: ***Ephemeris error at %.8fs\r\n', ...
                                    channel.PRN, channel.dataIndex/sampleFreq);
                            %-------------------------------------------------------------%
                            bits = -frameBuff(62) * frameBuff; %电平翻转
                            bits = dec2bin(bits>0)'; %±1数组转化为01字符串
                            fprintf(logID, ['%2d: ',bits(1:2)], channel.PRN);
                            for k=1:50 %将错误的电文输出，查找电文错误原因
                                fprintf(logID, ['%2d: ',bits((k-1)*30+2+(1:30))], channel.PRN);
                            end
                            %-------------------------------------------------------------%
                            frameBuffPoint = 0;
                            msgStage = 'H'; %重新寻找帧头
                            fprintf(logID, '%2d: Start find head at %.8fs\r\n', ...
                                    channel.PRN, channel.dataIndex/sampleFreq);
                        else
                            if ephemeris(2)~=ephemeris(3) %星历改变
                                fprintf(logID, '%2d: ***Ephemeris changes at %.8fs, IODC=%d, IODE=%d\r\n', ...
                                        channel.PRN, channel.dataIndex/sampleFreq, ephemeris(2), ephemeris(3));
                            else
                                channel.ephemeris = ephemeris; %更新星历
                                channel.state = 2; %更新状态（知道码发射时间，而且有星历）
                                fprintf(logID, '%2d: Ephemeris is parsed at %.8fs\r\n', ...
                                        channel.PRN, channel.dataIndex/sampleFreq);
                            end
                            frameBuff(1:2) = frameBuff(1501:1502); %将最后两个比特提前
                            frameBuffPoint = 2;
                        end
                    end
                otherwise
            end
        end
        
end

%% 更新下一数据块位置
trackDataTail = channel.trackDataHead + 1;
if trackDataTail>buffSize
    trackDataTail = 1;
end
blkSize = ceil((codeInt-remCodePhase)/codeNco*sampleFreq);
trackDataHead = trackDataTail + blkSize - 1;
if trackDataHead>buffSize
    trackDataHead = trackDataHead - buffSize;
end
channel.trackDataTail = trackDataTail;
channel.blkSize       = blkSize;
channel.trackDataHead = trackDataHead;

%% 更新通道信息
channel.trackStage     = trackStage;
channel.msgStage       = msgStage;
channel.cnt_t          = cnt_t;
channel.cnt_m          = cnt_m;
channel.code           = code;
channel.timeIntMs      = timeIntMs;
channel.carrNco        = carrNco;
channel.codeNco        = codeNco;
channel.carrAcc        = carrAcc;
channel.carrFreq       = carrFreq;
channel.codeFreq       = codeFreq;
channel.remCarrPhase   = remCarrPhase;
channel.remCodePhase   = remCodePhase;
channel.carrCirc       = carrCirc;
channel.I_P0           = I_P;
channel.Q_P0           = Q_P;
channel.FLL            = FLL;
channel.PLL            = PLL;
channel.DLL            = DLL;
channel.bitSyncTable   = bitSyncTable;
channel.bitBuff        = bitBuff;
channel.frameBuff      = frameBuff;
channel.frameBuffPoint = frameBuffPoint;
channel.CN0            = CN0;
channel.Px             = P;

%% 输出
I_Q = [I_P, I_E, I_L, Q_P, Q_E, Q_L];
disc = [codeError, codeSigma, carrError, carrSigma, freqError];

end