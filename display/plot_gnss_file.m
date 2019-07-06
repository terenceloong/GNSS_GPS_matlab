function plot_gnss_file(file_path)
% Plot first 0.1s data in specific file.

n = 4e5; %0.1s

fileID = fopen(file_path, 'r');
    data = fread(fileID, [2,n], 'int16'); %two row vector
    figure
    plot((1:n)/4e6, data(1,:))
    hold on
    plot((1:n)/4e6, data(2,:))
fclose(fileID);

end