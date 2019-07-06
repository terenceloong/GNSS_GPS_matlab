function print_log(filename, svList)
svList = svList'; %转成行向量，为了进行循环

fileID = fopen(filename, 'r');

for PRN=svList
    fprintf('PRN %d\n', PRN); %使用\r\n会多一个空行
    
    fseek(fileID, 0, 'bof'); %从文件头开始
    while ~feof(fileID) %一直读到文件尾
        tline = fgetl(fileID); %读一行
        [la, lb] = strtok(tline, ':'); %从:处分成两段
        ID = sscanf(la,'%2d'); %识别日志中的卫星编号
        if ID==PRN
            fprintf([lb(3:end),'\n']);
        end
    end
    
	disp(' ');
end

fclose(fileID);

end