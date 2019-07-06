function [a, xm] = mean_rec(a ,x1)
% 递推计算数据均值
% a.buff
% a.buffSize
% a.buffPoint，初值是0
% a.E0

k = a.buffSize; %缓存空间大小
n = a.buffPoint + 1; %指向当前缓存位置
x0 = a.buff(n); %当前缓存中的数据

a.E0 = a.E0 + (x1-x0)/k;

a.buff(n) = x1;
if n==k
    a.buffPoint = 0;
else
    a.buffPoint = n;
end

xm = a.E0;

end