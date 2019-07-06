function ch = GPS_L1_CA_channel_init(ch, acqResult, n, sampleFreq)
% 通道结构体初始化，执行完后通道被激活
% n表示已经经过了多少个采样点
% 需要预先赋PRN

code = GPS_L1_CA_generate(ch.PRN);

% ch.PRN 卫星号不变
ch.state = 1; %激活通道
ch.trackStage = 'F'; %频率牵引
ch.msgStage = 'I'; %空闲
ch.cnt_t = 0;
ch.cnt_m = 0;
ch.code = [code(end),code,code(1)]'; %列向量，为了求积分时用矢量相乘加速
ch.timeIntMs = 1;
ch.trackDataTail = sampleFreq*0.001 - acqResult(1) + 2;
ch.blkSize = sampleFreq*0.001;
ch.trackDataHead = ch.trackDataTail + ch.blkSize - 1;
ch.dataIndex = ch.trackDataTail + n;
ch.ts0 = NaN;
ch.carrNco = acqResult(2);
ch.codeNco = 1.023e6 + ch.carrNco/1540;
ch.carrAcc = 0;
ch.carrFreq = ch.carrNco;
ch.codeFreq = ch.codeNco;
ch.remCarrPhase = 0;
ch.remCodePhase = 0;
ch.carrCirc = 0;
ch.I_P0 = NaN;
ch.Q_P0 = NaN;

ch.FLL.K = 40;
ch.FLL.Int = ch.carrNco;

[K1, K2] = orderTwoLoopCoef(25, 0.707, 1);
ch.PLL.K1 = K1;
ch.PLL.K2 = K2;
ch.PLL.Int = 0;

[K1, K2] = orderTwoLoopCoef(2, 0.707, 1);
ch.DLL.K1 = K1;
ch.DLL.K2 = K2;
ch.DLL.Int = ch.codeNco;

ch.bitSyncTable = zeros(1,20);
ch.bitBuff = zeros(2,20); %第一行I_P，第二行Q_P
ch.frameBuff = zeros(1,1502);
ch.frameBuffPoint = 0;
ch.inverseFlag = -1; %-1表示不翻转，1表示翻转
% ch.ephemeris 星历不变

% 计算码鉴相器误差标准差结构体
ch.codeStd.buff = zeros(1,200);
ch.codeStd.buffSize = length(ch.codeStd.buff);
ch.codeStd.buffPoint = 0;
ch.codeStd.E0 = 0;
ch.codeStd.D0 = 0;

% 计算载波鉴相器误差标准差结构体
ch.carrStd.buff = zeros(1,200);
ch.carrStd.buffSize = length(ch.carrStd.buff);
ch.carrStd.buffPoint = 0;
ch.carrStd.E0 = 0;
ch.carrStd.D0 = 0;

% 计算NBP/WBP均值结构体
ch.NWmean.buff = zeros(1,50); %50个点求均值
ch.NWmean.buffSize = length(ch.NWmean.buff);
ch.NWmean.buffPoint = 0;
ch.NWmean.E0 = 0;
ch.CN0 = 0;

% 跟踪卡尔曼滤波器P阵
ch.Px = diag([0.02, 0.01, 5, 1].^2); %6m, 3.6deg, 5Hz, 1Hz/s (1sigma)

end