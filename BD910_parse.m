% Parse BD910 binary file.
clear;clc;

%% 1.Read file
dt = 100; %接收机输出时间间隔，ms

filename = './data/7_5/BD910_20190705_152853.dat';

fileID = fopen(filename, 'r');
stream = fread(fileID, 'uint8=>uint8');
fclose(fileID);
disp(['Read file: ',filename]);

%% 2.Delete the incomplete packet in the end
if stream(end)~=3
    disp('The last packet is incomplete!');
    while 1
        if (stream(end)==2)&&(stream(end-1)==3)
            stream(end)=[];
            break;
        else
            stream(end)=[];
        end
    end
end

%% 3.Count up the number of packets
n = length(stream);
pn = 0; %the number of packets
pp = 1;
while 1
    if stream(pp)==2 %head
        d = double(stream(pp+3)); %length
        if stream(pp+d+5)==3 %tail
            %----------------------------------------------%
            pn = pn+1;
            %----------------------------------------------%
            pp = pp+d+5; %pp points tail
        end
    end
    pp = pp+1; %pp points head
    if pp>n
        break;
    end
end
disp(['Packet account: ',num2str(pn)]);

%% 4.Checkout packet type
packet_type = zeros(pn,4);
k = 1;
pp = 1;
while 1
    if stream(pp)==2 %head
        d = double(stream(pp+3)); %length
        if stream(pp+d+5)==3 %tail
            %----------------------------------------------%
            packet_type(k,1) = double(stream(pp+2)); %packet type
            packet_type(k,2) = double(stream(pp+4)); %subtype
            packet_type(k,3) = pp;
            packet_type(k,4) = pp+d+5;
            k = k+1;
            %----------------------------------------------%
            pp = pp+d+5;
        end
    end
    pp = pp+1;
    if pp>n
        break;
    end
end

%% 5.Checkout whether lose packets
kn = 0; %the number of lost packets
for k=1:pn-1
    if packet_type(k+1,3)-packet_type(k,4)~=1
        disp(['Lose packet at packet ',num2str(k),'!']);
        kn = kn+1;
    end
end
if kn==0
    disp('No lost packet.');
end

%% 6.Parse position
disp('Parse position...');
fn = sum((packet_type(:,1)==87)&(packet_type(:,2)==1));
pos = zeros(fn,13);
k = 1;
pp = 1;
while 1
    if stream(pp)==2 %head
        d = double(stream(pp+3)); %length
        if stream(pp+d+5)==3 %tail
            %----------------------------------------------%
            if (stream(pp+2)==87)&&(stream(pp+4)==1)
                pos(k,:) = parse_position(stream(pp:pp+d+5));
                k = k+1;
            end
            %----------------------------------------------%
            pp = pp+d+5;
        end
    end
    pp = pp+1;
    if pp>n
        break;
    end
end
% fill lost data
k = 1;
while 1
    if pos(k,11)+dt~=pos(k+1,11)
        pos = [pos(1:k,:); NaN(1,13); pos(k+1:end,:)];
        pos(k+1,11) = pos(k,11)+dt;
        disp(['Fill data at pos ',num2str(k+1)]);
    end
    k = k+1;
    if k==length(pos)
        break;
    end
end

%% 7.Parse raw data
disp('Parse raw data...');
fn = sum((packet_type(:,1)==87)&(packet_type(:,2)==0));
info = zeros(fn,4);
raw = zeros(fn,32,8);
k = 1;
pp = 1;
while 1
    if stream(pp)==2 %head
        d = double(stream(pp+3)); %length
        if stream(pp+d+5)==3 %tail
            %----------------------------------------------%
            if (stream(pp+2)==87)&&(stream(pp+4)==0)
                [info(k,:), raw(k,:,:)] = parse_raw(stream(pp:pp+d+5));
                k = k+1;
            end
            %----------------------------------------------%
            pp = pp+d+5;
        end
    end
    pp = pp+1;
    if pp>n
        break;
    end
end
% fill lost data
k = 1;
while 1
    if info(k,2)==info(k+1,2)
        error(['Exception at raw ',num2str(k),'!']);
    end
    if info(k,2)+dt~=info(k+1,2)
        info = [info(1:k,:); NaN(1,4); info(k+1:end,:)];
        info(k+1,2) = info(k,2)+dt;
        raw = [raw(1:k,:,:); NaN(1,32,8); raw(k+1:end,:,:)];
        disp(['Fill data at raw ',num2str(k+1)]);
    end
    k = k+1;
    if k==length(info)
        break;
    end
end
flag1   = raw(:,:,1);
flag2   = raw(:,:,2);
ele     = raw(:,:,3);
azi     = raw(:,:,4);
snr     = raw(:,:,5);
rho     = raw(:,:,6);
phase   = raw(:,:,7);
doppler = raw(:,:,8);
clearvars raw

%% 8.Parse ephemeris
disp('Parse ephemeris...');
ion = zeros(0,9);
pp = 1;
while 1
    if stream(pp)==2 %head
        d = double(stream(pp+3)); %length
        if stream(pp+d+5)==3 %tail
            %----------------------------------------------%
            if (stream(pp+2)==87)&&(stream(pp+4)==0) %raw data
                t = typecast(swapbytes(typecast(stream(pp+(8:15)),'uint64')), 'double');
            end
            if (stream(pp+2)==85)&&(stream(pp+4)==3)
                eval('ion = [ion; [t,parse_ion(stream(pp:pp+d+5))]];');
            end
            if (stream(pp+2)==85)&&(stream(pp+4)==1)
                id = num2str(stream(pp+5));
                if ~exist(['ephemeris_',id], 'var')
                    eval(['ephemeris_',id,' = zeros(0,27);']);
                end
                eval(['ephemeris_',id,' = [ephemeris_',id,'; [t,parse_ephemeris(stream(pp:pp+d+5))]];']);
            end
            %----------------------------------------------%
            pp = pp+d+5;
        end
    end
    pp = pp+1;
    if pp>n
        break;
    end
end
disp('Finish!');

%%
clearvars -except packet_type pos info ...
    flag1 flag2 ele azi snr rho phase doppler ...
    ephemeris_* ion

% BD910_data.packet_type = packet_type;
% BD910_data.pos = pos;
% BD910_data.info = info;
% BD910_data.flag1 = flag1;
% BD910_data.flag2 = flag2;
% BD910_data.ele = ele;
% BD910_data.azi = azi;
% BD910_data.snr = snr;
% BD910_data.rho = rho;
% BD910_data.phase = phase;
% BD910_data.doppler = doppler;
% for k=1:32
%     ephemeris_name = ['ephemeris_',num2str(k)];
%     if exist(ephemeris_name, 'var')
%         eval(['BD910_data.',ephemeris_name,' = ',ephemeris_name,';'])
%     end
% end
% BD910_data.ion = ion;
% 
% clearvars -except BD910_data