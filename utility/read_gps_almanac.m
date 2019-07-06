function almanac = read_gps_almanac(filename)
%read GPS YUMA almanac file (.txt)

fileID = fopen(filename);
if(fileID == -1)
    almanac = [];
    disp('Can''t open the file!')
    return
end

temp = zeros(32,11);
almanac = zeros(32,12);
line = 1;
while ~feof(fileID)
    tline = fgetl(fileID);
    column = mod(line,15);
    if column>3 && column<15
        [~, remain] = strtok(tline, ':');
        C = textscan(remain, '%s %f');
        temp(row,column-3) = C{2};
    elseif column == 2
        [~, remain] = strtok(tline, ':');
        C = textscan(remain, '%s %f');
        row = C{2};
    end
    line = line+1;
end

temp(:,5) = temp(:,5).^2;
almanac(:,1) = (1:32)';          %ID
almanac(:,2) = temp(:,5);        %semi-major axis(m):a
almanac(:,3) = temp(:,1);        %eccentricity:e
almanac(:,4) = temp(:,3);        %orbital inclination(rad):i
almanac(:,5) = temp(:,6);        %ascending node right ascension(rad):Omega0
almanac(:,6) = temp(:,7);        %argument of perigee(rad):omega
almanac(:,7) = temp(:,8);        %mean anom(rad):M0
almanac(:,8) = temp(:,2);        %almanac reference time(s):toe
almanac(:,9) = temp(:,4);        %rate of right ascension(rad/s):OmegaDot
almanac(:,10) = temp(:,9);       %clock bias(s)
almanac(:,11) = temp(:,10);      %clock drift(s/s)
almanac(:,12) = temp(:,11);      %GPS week from 22 Aug 1999
n = 32;
i = 1;
while i <= n
    if almanac(i,2) == 0
        almanac(i,:) = [];
        n = n-1;
    else
        i = i+1;
    end
end

fclose(fileID);

end