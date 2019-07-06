function [gps_week, gps_second] = gps_time(c)
% Convert UTC/GMT+8 to GPS week and GPS second.
% clock = [year, month, date, hour, minute, second], UTC/GMT+8

% day = datenum(c(1),c(2),c(3)) - datenum(1999,8,22);
% day = datenum(c(1),c(2),c(3)) - datenum(2019,4,7);
if c(1)*10000+c(2)*100+c(3)>=20190407
    day = datenum(c(1),c(2),c(3)) - datenum(2019,4,7);
else
    day = datenum(c(1),c(2),c(3)) - datenum(1999,8,22);
end
gps_week = floor(day/7);
gps_second = (day-gps_week*7)*24*3600 + c(4)*3600 + c(5)*60 + floor(c(6));
gps_second = gps_second - 8*3600 + 18; %GPS time is 18 seconds faster than UTC.
if gps_second<0
    gps_second = gps_second+604800;
    gps_week = gps_week-1;
end

end