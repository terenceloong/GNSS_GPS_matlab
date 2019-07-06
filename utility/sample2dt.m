function dt = sample2dt(n, sampleFreq)
% 采样点数转化为时间增量（n必须大于0）

dt = [0,0,0]; %[s,ms,us]

% dt(1) = floor(n/sampleFreq);
% dt(2) = floor(rem(n,sampleFreq) * (1e3/sampleFreq));
% % (1e3/sampleFreq)表示一个采样点多少毫秒，rem(n,sampleFreq)表示不足1秒有多少个采样点
% dt(3) = rem(n,(sampleFreq/1e3)) * (1e6/sampleFreq);
% % (1e6/sampleFreq)表示一个采样点多少微秒，rem(n,(sampleFreq/1e3))表示不足1毫秒有多少个采样点

t = n / sampleFreq;
dt(1) = floor(t); %整秒部分
t = mod(t,1) * 1000;
dt(2) = floor(t); %毫秒部分
dt(3) = mod(t,1) * 1000; %微秒部分

end