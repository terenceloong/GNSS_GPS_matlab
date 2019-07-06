function trackResult = trackResult_struct(m)
% 跟踪结果结构体

trackResult.PRN = 0;
trackResult.n = 1; %指向当前存储的行号

trackResult.dataIndex     = zeros(m,1); %码周期开始采样点在原始数据文件中的位置
trackResult.ts0           = zeros(m,1); %码周期理论发射时间，ms
trackResult.remCodePhase  = zeros(m,1); %码周期开始采样点的码相位，码片
trackResult.codeFreq      = zeros(m,1); %码频率
trackResult.remCarrPhase  = zeros(m,1); %码周期开始采样点的载波相位，周
trackResult.carrFreq      = zeros(m,1); %载波频率
trackResult.I_Q           = zeros(m,6); %[I_P,I_E,I_L,Q_P,Q_E,Q_L]
trackResult.disc          = zeros(m,5); %[codeError,std, carrError,std, freqError]，鉴相器
trackResult.bitStartFlag  = zeros(m,1); %比特开始标志
trackResult.CN0           = zeros(m,1); %载噪比
trackResult.carrAcc       = zeros(m,1); %载波加速度
trackResult.Px            = zeros(m,4); %跟踪滤波器的P阵

end