function ephemeris = parse_ephemeris(stream)

ephemeris = zeros(1,26);

stream(1) = [];
ephemeris(1)  = double(typecast(swapbytes(typecast(stream(6:7),'uint16')), 'uint16')); %week
ephemeris(2)  = double(typecast(swapbytes(typecast(stream(8:9),'uint16')), 'uint16')); %IODC
ephemeris(3)  = double(stream(11)); %IODE
ephemeris(4)  = double(typecast(swapbytes(typecast(stream(12:15),'uint32')), 'uint32')); %TOW
ephemeris(5)  = double(typecast(swapbytes(typecast(stream(16:19),'uint32')), 'uint32')); %toc
ephemeris(6)  = double(typecast(swapbytes(typecast(stream(20:23),'uint32')), 'uint32')); %toe
ephemeris(7)  = typecast(swapbytes(typecast(stream(24:31),'uint64')), 'double'); %tGD
ephemeris(8)  = typecast(swapbytes(typecast(stream(32:39),'uint64')), 'double'); %af2
ephemeris(9)  = typecast(swapbytes(typecast(stream(40:47),'uint64')), 'double'); %af1
ephemeris(10) = typecast(swapbytes(typecast(stream(48:55),'uint64')), 'double'); %af0
ephemeris(11) = typecast(swapbytes(typecast(stream(56:63),'uint64')), 'double'); %Crs
ephemeris(12) = typecast(swapbytes(typecast(stream(64:71),'uint64')), 'double') * pi; %dn
ephemeris(13) = typecast(swapbytes(typecast(stream(72:79),'uint64')), 'double') * pi; %M0
ephemeris(14) = typecast(swapbytes(typecast(stream(80:87),'uint64')), 'double') * pi; %Cuc
ephemeris(15) = typecast(swapbytes(typecast(stream(88:95),'uint64')), 'double'); %e
ephemeris(16) = typecast(swapbytes(typecast(stream(96:103),'uint64')), 'double') * pi; %Cus
ephemeris(17) = typecast(swapbytes(typecast(stream(104:111),'uint64')), 'double'); %sqa
ephemeris(18) = typecast(swapbytes(typecast(stream(112:119),'uint64')), 'double') * pi; %Cic
ephemeris(19) = typecast(swapbytes(typecast(stream(120:127),'uint64')), 'double') * pi; %Omega0
ephemeris(20) = typecast(swapbytes(typecast(stream(128:135),'uint64')), 'double') * pi; %Cis
ephemeris(21) = typecast(swapbytes(typecast(stream(136:143),'uint64')), 'double') * pi; %i0
ephemeris(22) = typecast(swapbytes(typecast(stream(144:151),'uint64')), 'double'); %Crc
ephemeris(23) = typecast(swapbytes(typecast(stream(152:159),'uint64')), 'double') * pi; %omega
ephemeris(24) = typecast(swapbytes(typecast(stream(160:167),'uint64')), 'double') * pi; %Omega_dot
ephemeris(25) = typecast(swapbytes(typecast(stream(168:175),'uint64')), 'double') * pi; %i_dot
ephemeris(26) = double(typecast(swapbytes(typecast(stream(176:179),'uint32')), 'uint32')); %flags

end