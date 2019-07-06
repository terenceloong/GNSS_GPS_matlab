function [info, raw] = parse_raw(stream)

info = zeros(1,4);
raw = NaN(32,8);

stream(1) = [];
info(1) = double(stream(6)); %reply number
info(2) = typecast(swapbytes(typecast(stream(8:15),'uint64')), 'double'); %receive time, ms
info(3) = typecast(swapbytes(typecast(stream(16:23),'uint64')), 'double'); %clock offset, ms
info(4) = double(stream(24)); %# of SVs in record
d = double(stream(3)); %length
m = 0;
while m<d-21 %25+m<d+5-1
    id = stream(25+m);
    raw(id,1) = double(stream(26+m)); %flag1
    raw(id,2) = double(stream(27+m)); %flag2
    raw(id,3) = double(stream(28+m)); %elevation angle
    raw(id,4) = double(typecast(swapbytes(typecast(stream((29:30)+m),'uint16')), 'int16')); %azimuth
    raw(id,5) = double(stream(31+m))/4; %SNR
    raw(id,6) = typecast(swapbytes(typecast(stream((32:39)+m),'uint64')), 'double'); %pseudo range
    raw(id,7) = typecast(swapbytes(typecast(stream((40:47)+m),'uint64')), 'double'); %phase
    raw(id,8) = double(typecast(swapbytes(typecast(stream((48:51)+m),'uint32')), 'single')); %doppler
    m = m+27;
end

end