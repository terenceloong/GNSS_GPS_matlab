function [acqResult, peakRatio] = GPS_L1_CA_acq_one(PRN, data)
% 捕获编号为PRN的卫星。如果未捕获到，acqResult为空
% data需要有连续两个捕获长度的数据

acqResult = []; %存捕获结果，第一列码相位，第二列载波频率

%%
N = length(data) / 2; %采样点数
fs = 4e6; %采样频率，Hz
fc = 1.023e6; %码频率，Hz

carrFreq = -6e3:(fs/N/2):6e3; %频率搜索范围
M = length(carrFreq);

acqThreshold = 1.4; %搜索阈值，最大峰比第二大峰大多少倍

%%
% 取连续两个复基带数据，因为可能存在导航电文的翻转导致相关峰变小
baseband1 = data(1,  1:N)   + data(2,  1:N)  *1i;
baseband2 = data(1,N+1:end) + data(2,N+1:end)*1i;

result1 = zeros(M,4000); %存搜索结果，行是载波频率，列是码相位
result2 = zeros(M,4000);
resultCorr1 = zeros(M,2); %第一列存每个搜索频率相关最大值，第二列存最大值对应的索引
resultCorr2 = zeros(M,2);

% [Xg,Yg] = meshgrid(1:4000,carrFreq); %画图格子-------------

%% 捕获算法
CAcode = GPS_L1_CA_generate(PRN); %C/A码序列
code = CAcode(mod(floor((0:N-1)*fc/fs),1023) + 1); %C/A码采样
CODE = fft(code); %C/A码FFT

%----搜索
for k=1:M
    carrier = exp(-2*pi * carrFreq(k) * (0:N-1)/fs * 1i); %本地复载波，负频率

    x = baseband1 .* carrier;
    X = fft(x);
    Y = conj(X).*CODE;
    y = abs(ifft(Y));
    result1(k,:) = y(1:4000); %只取前4000个数，后面的都是重复的
    [resultCorr1(k,1), resultCorr1(k,2)] = max(result1(k,:)); %寻找一次相关的最大值及其索引

    x = baseband2 .* carrier;
    X = fft(x);
    Y = conj(X).*CODE;
    y = abs(ifft(Y));
    result2(k,:) = y(1:4000);
    [resultCorr2(k,1), resultCorr2(k,2)] = max(result2(k,:));
end

%----选取值大的那组数据
if max(resultCorr1(:,1))>max(resultCorr2(:,1))
    corrValue = resultCorr1(:,1);
    corrIndex = resultCorr1(:,2);
%     result = result1; %用来画图-------------
else
    corrValue = resultCorr2(:,1);
    corrIndex = resultCorr2(:,2);
%     result = result2; %用来画图-------------
end

%----寻找相关峰
[peakSize, index] = max(corrValue); %最大峰
corrValue(mod(index+(-3:3)-1,M)+1) = 0; %排除掉最大相关峰周围的点
secondPeakSize = max(corrValue); %第二大峰

%----捕获到信号
peakRatio = peakSize / secondPeakSize; %最高峰与第二大峰的比值
if peakRatio>acqThreshold
    % 画图-------------
%     figure
%     surf(Xg,Yg,result)
%     title(['PRN = ',num2str(PRN)])

    acqResult(1) = corrIndex(index); %码相位
    acqResult(2) = carrFreq(index); %载波频率
end

end