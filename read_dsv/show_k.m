%
% Plot some k-space from simulation files
%


ky = Read_dsv('C:/Temp/DspData_M0Y.dsv');
kx = Read_dsv('C:/Temp/DspData_M0X.dsv');

ADC = Read_dsv('C:/Temp/DspData_ADC.dsv');

plot(kx.timecourse(points),ky.timecourse(points))