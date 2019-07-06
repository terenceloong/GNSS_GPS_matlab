function [filename, gps_week, gps_second] = download_gps_almanac(c)
% 下载GPS历书，提前建好almanac文件夹
% c为本地时间，[year, month, day, hour, minute, second]

[gps_week, gps_second] = gps_time(c);

if gps_second<61440
    w = sprintf('%04d',gps_week);
    s = '061440';
elseif gps_second<147456
    w = sprintf('%04d',gps_week);
    s = '147456';
elseif gps_second<233472
    w = sprintf('%04d',gps_week);
    s = '233472';
elseif gps_second<319488
    w = sprintf('%04d',gps_week);
    s = '319488';
elseif gps_second<405504
    w = sprintf('%04d',gps_week);
    s = '405504';
elseif gps_second<589824
    w = sprintf('%04d',gps_week);
    s = '589824';
else
    w = sprintf('%04d',gps_week+1);
    s = '061440';
end
filename = ['./almanac/',w,'_',s,'.txt'];

if ~exist(['./almanac/',w,'_',s,'.txt'], 'file')
    websave(filename, ['http://celestrak.com/GPS/almanac/Yuma/',num2str(c(1)),'/almanac.yuma.week',w,'.',s,'.txt']);
end

end