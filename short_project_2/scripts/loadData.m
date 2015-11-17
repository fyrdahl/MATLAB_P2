load([savingdir,'/MAT-files/images/',phantomName,'_list',num2str(offsetListNum),'TEimages.mat'],'TEimages')
load([savingdir,'/MAT-files/images/',phantomName,'_list',num2str(offsetListNum),'TIimages.mat'],'TIimages')
load([savingdir,'/MAT-files/images/',phantomName,'_list',num2str(offsetListNum),'fullFPimages.mat'],'fullFPimages')
load([savingdir,'/MAT-files/images/',phantomName,'_list',num2str(offsetListNum),'TE.mat'],'TE')
load([savingdir,'/MAT-files/images/',phantomName,'_list',num2str(offsetListNum),'TI.mat'],'TI')
FPimages(:,:,:) = squeeze(fullFPimages(:,:,sliceNumber,:));