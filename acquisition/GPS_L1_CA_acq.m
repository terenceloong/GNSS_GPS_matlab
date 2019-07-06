function acqResults = GPS_L1_CA_acq(file_path, sample_offset, N)
% GPS信号捕获，32颗卫星全搜索，连续搜索2段数据，结果存在变量acqResults中
% sample_offset：抛弃前多少个采样点处开始处理
% N：FFT点数

%%
Ns = N; %采样点数
fs = 4e6; %采样频率，Hz
fc = 1.023e6; %码频率，Hz

carrFreq = -6e3:(fs/Ns/2):6e3; %频率搜索范围
M = length(carrFreq);

acqThreshold = 1.4; %搜索阈值，最大峰比第二大峰大多少倍

%%
% 取连续两个复基带数据，因为可能存在导航电文的翻转导致相关峰变小
fileID = fopen(file_path, 'r');
    fseek(fileID, round(sample_offset*4), 'bof');
    if int32(ftell(fileID))~=int32(sample_offset*4)
        error('Sample offset error!');
    end
    baseband1 = double(fread(fileID, [2,Ns], 'int16'));
    baseband1 = baseband1(1,:) + baseband1(2,:)*1i; %行向量
    baseband2 = double(fread(fileID, [2,Ns], 'int16'));
    baseband2 = baseband2(1,:) + baseband2(2,:)*1i; %行向量
fclose(fileID);

result1 = zeros(M,4000); %存搜索结果，行是载波频率，列是码相位
result2 = zeros(M,4000);
resultCorr1 = zeros(M,2); %第一列存每个搜索频率相关最大值，第二列存最大值对应的索引
resultCorr2 = zeros(M,2);

[Xg,Yg] = meshgrid(1:4000,carrFreq); %画图格子

% 存捕获结果，第一列码相位，第二列载波频率
acqResults = NaN(32,2);

%% 捕获算法
for PRN=1:32
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
        result = result1; %用来画图
    else
        corrValue = resultCorr2(:,1);
        corrIndex = resultCorr2(:,2);
        result = result2;
    end
    
    %----寻找相关峰
    [peakSize, index] = max(corrValue); %最大峰
    corrValue(mod(index+(-3:3)-1,M)+1) = 0; %排除掉最大相关峰周围的点
    secondPeakSize = max(corrValue); %第二大峰
    
    %----捕获到信号
    if (peakSize/secondPeakSize)>acqThreshold
        % 画图
        figure
        surf(Xg,Yg,result)
        title(['PRN = ',num2str(PRN)])
        
        %存储捕获结果
        acqResults(PRN,1) = corrIndex(index); %码相位
        acqResults(PRN,2) = carrFreq(index); %载波频率
    end
end

end