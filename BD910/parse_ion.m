function ion = parse_ion(stream)

ion = zeros(1,8);

stream(1) = [];
ion(1) = typecast(swapbytes(typecast(stream(6 :13),'uint64')), 'double'); %alpha1
ion(2) = typecast(swapbytes(typecast(stream(14:21),'uint64')), 'double'); %alpha2
ion(3) = typecast(swapbytes(typecast(stream(22:29),'uint64')), 'double'); %alpha3
ion(4) = typecast(swapbytes(typecast(stream(30:37),'uint64')), 'double'); %alpha4
ion(5) = typecast(swapbytes(typecast(stream(38:45),'uint64')), 'double'); %beta1
ion(6) = typecast(swapbytes(typecast(stream(46:53),'uint64')), 'double'); %beta2
ion(7) = typecast(swapbytes(typecast(stream(54:61),'uint64')), 'double'); %beta3
ion(8) = typecast(swapbytes(typecast(stream(62:69),'uint64')), 'double'); %beta4

end