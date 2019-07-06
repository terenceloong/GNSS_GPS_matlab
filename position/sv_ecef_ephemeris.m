function [sv, E] = sv_ecef_ephemeris(ephemeris, t)
% 给定GPS时间t计算卫星在地心固连坐标系下的位置和速度
% 输出偏近点角E，用来后续算相对论校正
% sv = [x,y,z, vx,vy,vz]

miu = 3.986005e14;
w = 7.2921151467e-5;

sqa = ephemeris(17);
toe = ephemeris(6);
dn = ephemeris(12);
M0 = ephemeris(13);
e = ephemeris(15);
omega = ephemeris(23);
Cus = ephemeris(16);
Cuc = ephemeris(14);
Crs = ephemeris(11);
Crc = ephemeris(22);
Cis = ephemeris(20);
Cic = ephemeris(18);
i0 = ephemeris(21);
i_dot = ephemeris(25);
Omega0 = ephemeris(19);
Omega_dot = ephemeris(24);

% 参见《ICD-GPS》
a = sqa^2;
n0 = sqrt(miu/a^3);
dt = t - toe;
if dt>302400
    dt = dt-604800;
elseif dt<-302400
    dt = dt+604800;
end
n = n0 + dn;
M = mod(M0+n*dt, 2*pi); %平近点角，0~2*pi
E = kepler(M, e); %偏近点角
v = 2 * mod(atan(sqrt((1+e)/(1-e))*tan(E/2)), pi); %真近点角，0~2*pi
Phi = v + omega;
du = Cus*sin(2*Phi) + Cuc*cos(2*Phi);
dr = Crs*sin(2*Phi) + Crc*cos(2*Phi);
di = Cis*sin(2*Phi) + Cic*cos(2*Phi);
u = Phi + du;
r = a*(1-e*cos(E)) + dr;
i = i0 + di + i_dot*dt;
xp = r*cos(u);
yp = r*sin(u);
Omega = Omega0 + (Omega_dot-w)*dt - w*toe;
x = xp*cos(Omega) - yp*cos(i)*sin(Omega);
y = xp*sin(Omega) + yp*cos(i)*cos(Omega);
z = yp*sin(i);

% 参见《北斗/GPS双模软件接收机原理与实现技术》283页
d_E = n/(1-e*cos(E));
d_Phi = sqrt(1-e^2)*d_E/(1-e*cos(E));
d_r = a*e*sin(E)*d_E + 2*(Crs*cos(2*Phi)-Crc*sin(2*Phi))*d_Phi;
d_u = d_Phi + 2*(Cus*cos(2*Phi)-Cuc*sin(2*Phi))*d_Phi;
d_Omega = Omega_dot-w;
d_i = i_dot + 2*(Cis*cos(2*Phi)-Cic*sin(2*Phi))*d_Phi;
d_xp = d_r*cos(u) - r*sin(u)*d_u;
d_yp = d_r*sin(u) + r*cos(u)*d_u;
vx = d_xp*cos(Omega) - d_yp*cos(i)*sin(Omega) + yp*sin(i)*sin(Omega)*d_i - y*d_Omega;
vy = d_xp*sin(Omega) + d_yp*cos(i)*cos(Omega) - yp*sin(i)*cos(Omega)*d_i + x*d_Omega;
vz = d_yp*sin(i) + yp*cos(i)*d_i;

sv = zeros(1,6);
sv(1:3) = [x,y,z]; %position (ecef)
sv(4:6) = [vx,vy,vz]; %velocity (ecef)

end

function E = kepler(M, e)
    E = M;
    Ei = E - (E-e*sin(E)-M)/(1-e*cos(E));
    while abs(Ei-E) > 1e-12
        E = Ei;
        Ei = E - (E-e*sin(E)-M)/(1-e*cos(E));
    end
end