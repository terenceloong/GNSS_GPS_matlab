function calib_baseline()
% 标定基线/测姿态

%% 输入数据
svList = evalin('base', 'svList');
svN = evalin('base', 'svN');
output_pos = evalin('base', 'output_pos');
output_sv = evalin('base', 'output_sv');
output_dphase = evalin('base', 'output_dphase');

%% 参数
bl = 1.32; %大致基线长度
br = 0.02; %基线长度范围
n = size(output_pos,1); %计算多少个点
lamda = 299792458 / 1575.42e6; %波长

circ_limit = 1000;
circ_half = circ_limit/2;

mode = 0; %0为标定基线，1为测姿态

%% 结果变量
BLs = NaN(n,3); %基线测量结果，第一列航向角，第二列俯仰角，第三列基线长度
pdb = NaN(n,svN); %根据基线算的相位差（理论相位差）
pdm = NaN(n,svN); %实测的相位差加整周（实测相位差）

%% 计算基线矢量（全程整周模糊度搜索）
if mode==0
    for k=1:n
        %----接收机位置
        p0 = output_pos(k,1:3); %纬经高
        Cen = dcmecef2ned(p0(1), p0(2));
        rp = lla2ecef(p0); %ecef坐标，行向量
        %----地理系下卫星指向天线的单位矢量
        rs = output_sv(:,1:3,k); %卫星坐标
        rsp = ones(svN,1)*rp - rs;
        rho = sum(rsp.*rsp, 2).^0.5; %列向量
        rspu = rsp ./ (rho*[1,1,1]);
        A = rspu * Cen';
        %----相位差
        p = output_dphase(k,:)'; %列向量，通道输出的相位差
        p = mod(p,1); %取小数部分
        %----排除高度角太低的卫星
        % p(asind(A(:,3))<20) = NaN;
        %----删除没有相位差的行
        Ac = A(~isnan(p),:);
        pc = p(~isnan(p));
        %----计算基线矢量
        if length(pc)>=5
            Rx = IAR_nonbaseline(Ac, pc, lamda, bl+[-br,br]);
            L = norm(Rx); %基线长度
            psi = atan2d(Rx(2),Rx(1)); %基线航向角
            theta = -asind(Rx(3)/L); %基线俯仰角
            % 存储结果
            BLs(k,:) = [psi,theta,L];
            pdb(k,:) = (A*Rx / lamda)';
            pdm(k,:) = p' + round(pdb(k,:)-p');
        end
    end
end
    
%% 计算基线矢量（记录整周模糊度）
if mode==1
    N = NaN(svN,1); %所有通道的相位整周误差，通道输出的相位差减去这个值为真实相位差
    for k=1:n
        %----接收机位置
        p0 = output_pos(k,1:3); %纬经高
        Cen = dcmecef2ned(p0(1), p0(2));
        rp = lla2ecef(p0); %ecef坐标，行向量
        %----地理系下卫星指向天线的单位矢量
        rs = output_sv(:,1:3,k); %卫星坐标
        rsp = ones(svN,1)*rp - rs;
        rho = sum(rsp.*rsp, 2).^0.5; %列向量
        rspu = rsp ./ (rho*[1,1,1]);
        A = rspu * Cen';
        %----相位差
        p0 = output_dphase(k,:)'; %列向量，通道输出的相位差
        %----计算基线矢量
        if sum(~isnan(p0-N))<4 %完整的相位差数量小于4，不能定姿
            if sum(~isnan(p0))>=5 %测出的相位差数量大于等于5，可以进行整周模糊度搜索
                p = mod(p0,1); %取小数部分
                Ac = A(~isnan(p),:); %删除没有相位差的行
                pc = p(~isnan(p));
                Rx = IAR_nonbaseline(Ac, pc, lamda, bl+[-br,br]);
                N = round(p0-A*Rx/lamda); %计算所有通道的相位整周误差
            else
                continue
            end
        else %完整的相位差数量大于等于4，可以直接定姿
            p = mod(p0-N + circ_half, circ_limit) - circ_half; %修正后的相位差
            Ac = A(~isnan(p),:); %删除没有相位差的行
            pc = p(~isnan(p));
            ele = asind(Ac(:,3));
            [~,i1] = max(ele);
            ele(i1) = [];
            W = diag(sind(ele))^2; %权值
            Ac = Ac - ones(size(Ac,1),1)*Ac(i1,:);
            Ac(i1,:) = [];
            pc = pc - pc(i1);
            pc(i1) = [];
            Rx = (Ac'*W*Ac) \ (Ac'*W*pc*lamda); %加权最小二乘
            N = round(p0-A*Rx/lamda); %计算所有通道的相位整周误差，新捕获通道对应的值会被直接计算，中断通道的值会被销毁
        end
        %----存储结果
        L = norm(Rx); %基线长度
        psi = atan2d(Rx(2),Rx(1)); %基线航向角
        theta = -asind(Rx(3)/L); %基线俯仰角
        BLs(k,:) = [psi,theta,L];
        pdb(k,:) = (A*Rx / lamda)';
        pdm(k,:) = (mod(p0-N + circ_half, circ_limit) - circ_half)';
    end
end

%% 输出数据
assignin('base', 'BLs', BLs)
assignin('base', 'pdb', pdb)
assignin('base', 'pdm', pdm)

%% 画基线标定结果
figure
subplot(3,1,1)
plot(BLs(:,1))
grid on
title('航向角')
subplot(3,1,2)
plot(BLs(:,2))
grid on
title('俯仰角')
subplot(3,1,3)
plot(BLs(:,3))
grid on
title('基线长度')
if exist('br','var')
    set(gca,'Ylim',[bl-br,bl+br])
end

%% 画相位差
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
                  0,     0,     0];

figure
hold on
grid on
legend_str = [];
for k=1:svN
    if sum(~isnan(pdm(:,k)))~=0
        plot(pdm(:,k), 'Color',colorTable(k,:)) %实线，实测相位差
        legend_str = [legend_str; string(num2str(svList(k)))];
    end
end
for k=1:svN
    if sum(~isnan(pdm(:,k)))~=0
        plot(pdb(:,k), 'Color',colorTable(k,:), 'LineStyle','--') %虚线，基线算的相位差
    end
end
legend(legend_str)
title('实测相位差与计算相位差')

%% 画实测相位差与计算相位差之差（最好在一条直线上）
% 两天线线路不同造成的
pdd = pdm - pdb;
figure
hold on
grid on
for k=1:svN
    if sum(~isnan(pdd(:,k)))~=0
        plot(pdd(:,k), 'Color',colorTable(k,:))
    end
end
legend(legend_str)
set(gca,'Ylim',[-0.5,0.5])
title('实测相位差 - 计算相位差')

end