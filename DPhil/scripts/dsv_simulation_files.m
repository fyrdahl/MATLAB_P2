%% Example of reading in POET simulation data
addpath(genpath('~/Documents/MATLAB/read_dsv'));
% x-axis gradient
dsv_filename = '~/Documents/DPhil/dsv_files';
dsv = Read_dsv(dsv_filename);

output_data = dsv2timecourse(dsv);