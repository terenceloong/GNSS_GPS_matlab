function intNumber = twosComp2dec(binaryNumber)
% 转化二进制补码字符串为十进制数字

%--- Convert from binary form to a decimal number -------------------------
intNumber = bin2dec(binaryNumber);

%--- If the number was negative, then correct the result ------------------
if binaryNumber(1) == '1'
    intNumber = intNumber - 2^length(binaryNumber);
end

end