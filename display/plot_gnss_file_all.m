% 浏览文件中的所有数据

file_path = '.\data\data_20190704_210904_ch1.dat';

n = 4e5; %0.1s

fileID = fopen(file_path, 'r');

k = 0;
figure
while 1
    data = fread(fileID, [2,n], 'int16'); %two row vector
    if size(data,2)~=n
        break
    end
    plot((1:n)/4e6+0.1*k, data(1,:))
    hold on
    plot((1:n)/4e6+0.1*k, data(2,:))
    hold off
    set(gca, 'YLim',[-600,600])
    drawnow
    k = k+1;
end
    
fclose(fileID);