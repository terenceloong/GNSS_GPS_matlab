function pos = parse_position(stream)

pos = zeros(1,13);

stream(1) = [];
pos(1) = double(stream(6)); %reply number
pos(2) = typecast(swapbytes(typecast(stream(8 :15),'uint64')), 'double') * 180; %latitude, deg
pos(3) = typecast(swapbytes(typecast(stream(16:23),'uint64')), 'double') * 180; %longitude, deg
pos(4) = typecast(swapbytes(typecast(stream(24:31),'uint64')), 'double'); %altitude, m
pos(5) = typecast(swapbytes(typecast(stream(32:39),'uint64')), 'double'); %clock offset, m
pos(6) = typecast(swapbytes(typecast(stream(40:47),'uint64')), 'double'); %frequency offset, Hz
pos(7) = typecast(swapbytes(typecast(stream(48:55),'uint64')), 'double'); %PDOP
pos(8) = typecast(swapbytes(typecast(stream(56:63),'uint64')), 'double'); %latitude rate, rad/s
pos(9) = typecast(swapbytes(typecast(stream(64:71),'uint64')), 'double'); %longitude rate, rad/s
pos(10) = typecast(swapbytes(typecast(stream(72:79),'uint64')), 'double'); %altitude rade, rad/s
pos(11) = double(typecast(swapbytes(typecast(stream(80:83),'uint32')), 'uint32')); %GPS msec of week, ms
pos(12) = double(stream(84)); %flag
pos(13) = double(stream(85)); %# of SVs, number of satellites used to compute position solution

end