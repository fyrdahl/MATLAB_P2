function [FPimages] = load_FP_Data(phantomName, offsetListNum, sliceNumber, savingdir)
%Function name: loadData
%Description: Loads fingerprinting MR data from a specified directory.

load([savingdir,'/MAT-files/images/',phantomName,'_list',num2str(offsetListNum),'fullFPimages.mat'],'fullFPimages');

FPimages(:,:,:) = squeeze(fullFPimages(:,:,sliceNumber,:));

end