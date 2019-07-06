%% 单天线
clear; clc;
sample_offset = 0*4e6;
acqResults = GPS_L1_CA_acq('.\data\7_2\data_20190702_111609_ch1.dat', sample_offset, 8000);

%% 双天线
clear; clc;
sample_offset = 0*4e6;
acqResults_A = GPS_L1_CA_acq('.\data\7_2\data_20190702_111609_ch1.dat', sample_offset, 8000);
acqResults_B = GPS_L1_CA_acq('.\data\7_2\data_20190702_111609_ch2.dat', sample_offset, 8000);