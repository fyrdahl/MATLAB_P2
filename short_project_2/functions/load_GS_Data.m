function [TEimages, TIimages, TE, TI] = load_GS_Data(phantomName, offsetListNum, savingdir)

%[TEimages, TIimages, TE, TI] = loadGSData(phantomName, offsetListNum, sliceNumber, savingdir)
%
%Description: Loads gold standard MR data for measuring T1 and T2. TI
%images acquired with inversion recovery and TE images acquired with spin
%echo.
%
%Author: Jack Allen <jack.allen@jesus.ox.ac.uk>

load([savingdir,'/MAT-files/images/',phantomName,'_list',num2str(offsetListNum),'TEimages.mat'],'TEimages')
load([savingdir,'/MAT-files/images/',phantomName,'_list',num2str(offsetListNum),'TIimages.mat'],'TIimages')
load([savingdir,'/MAT-files/images/',phantomName,'_list',num2str(offsetListNum),'TE.mat'],'TE')
load([savingdir,'/MAT-files/images/',phantomName,'_list',num2str(offsetListNum),'TI.mat'],'TI')

end