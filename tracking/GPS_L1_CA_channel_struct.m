function ch = GPS_L1_CA_channel_struct()
% 声明通道结构体所有场

ch.PRN              = []; %卫星编号
ch.state            = []; %通道状态（数字）
ch.trackStage       = []; %跟踪阶段（字符）
ch.msgStage         = []; %电文解析阶段（字符）
ch.cnt_t            = []; %跟踪时用的计数器
ch.cnt_m            = []; %电文解析时用的计数器
ch.code             = []; %伪码
ch.timeIntMs        = []; %积分时间，ms
ch.trackDataTail    = []; %跟踪开始点在数据缓存中的位置
ch.blkSize          = []; %跟踪数据段采样点个数
ch.trackDataHead    = []; %跟踪结束点在数据缓存中的位置
ch.dataIndex        = []; %跟踪开始点在文件中的位置
ch.ts0              = []; %跟踪开始点所在码周期的理论发射时间，ms
ch.carrNco          = []; %载波发生器频率
ch.codeNco          = []; %码发生器频率
ch.carrAcc          = []; %载波频率加速度
ch.carrFreq         = []; %载波频率测量
ch.codeFreq         = []; %码频率测量
ch.remCarrPhase     = []; %跟踪开始点的载波相位
ch.remCodePhase     = []; %跟踪开始点的码相位
ch.carrCirc         = []; %记录载波经过的整周数，0~999
ch.I_P0             = []; %上次跟踪的I_P
ch.Q_P0             = []; %上次跟踪的Q_P
ch.FLL              = []; %锁频环（结构体）
ch.PLL              = []; %锁相环（结构体）
ch.DLL              = []; %延迟锁定环（结构体）
ch.bitSyncTable     = []; %比特同步统计表
ch.bitBuff          = []; %比特缓存
ch.frameBuff        = []; %帧缓存
ch.frameBuffPoint   = []; %帧缓存指针
ch.inverseFlag      = []; %相位翻转标志
ch.ephemeris        = []; %星历
ch.codeStd          = []; %计算码鉴相器误差标准差结构体
ch.carrStd          = []; %计算载波鉴相器误差标准差结构体
ch.NWmean           = []; %计算NBP/WBP均值结构体
ch.CN0              = []; %载噪比
ch.Px               = []; %卡尔曼滤波跟踪方差阵

end