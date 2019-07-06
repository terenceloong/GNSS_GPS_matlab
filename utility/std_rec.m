function [a, sigma] = std_rec(a ,x1)
% 递推计算数据标准差
% a.buff
% a.buffSize
% a.buffPoint，初值是0
% a.E0
% a.D0

k = a.buffSize; %缓存空间大小
n = a.buffPoint + 1; %指向当前缓存位置
x0 = a.buff(n); %当前缓存中的数据

E0 = a.E0;
D0 = a.D0;

E1 = E0 + (x1-x0)/k;
D1 = D0 + ((x1-E1)^2 - (x0-E0)^2 - 2*(E1-E0)*(E0*k-x0) + (k-1)*(E1^2-E0^2))/k;

a.E0 = E1;
a.D0 = D1;
a.buff(n) = x1;
if n==k
    a.buffPoint = 0;
else
    a.buffPoint = n;
end

sigma = sqrt(D1);

end