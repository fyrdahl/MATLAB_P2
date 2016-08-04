
function [testT1map, testT2map] = make_synthetic_phantom(phantom, varargin )
%% Make Test T1 and T2 maps
switch phantom
    case 'brainTemplate'
        load ('~/Documents/MATLAB/FMRIB_graduate_course/3_Contrast_Manipulation/segBrain.mat')
        testT1map = zeros(size(segBrain,1), size(segBrain,2));
        testT2map = zeros(size(segBrain,1), size(segBrain,2));
        
        [T1, T2] = get_relaxation_times(3,'wm');
        testT1map(segBrain == 1) = T1;
        testT2map(segBrain == 1) = T2;
        [T1, T2] = get_relaxation_times(3,'gm');
        testT1map(segBrain == 2) = T1;
        testT2map(segBrain == 2) = T2;
        %CSF
        %Ma MRM2013
        testT1map(segBrain == 3) = 4880;
        testT2map(segBrain == 3) = 550;
        
    case 'sphere'
        T1 = 282.3;
        T2 = 214.8;
        dims = [varargin{1} varargin{2}];
        testT1map = T1*ones(dims(1), dims(2));
        testT2map = T2*ones(dims(1), dims(2));
        
        if dims(1)>=20 && dims(2)>=20
            testT1map(5:10,5:10) = 230;
            testT1map(5:10,10:15) = 260;
            testT1map(10:15,5:10) = 280;
            testT1map(10:15,10:15) = 300;
            
            testT2map(5:10,5:10) = 210;
            testT2map(5:10,10:15) = 220;
            testT2map(10:15,5:10) = 240;
            testT2map(10:15,10:15) = 280;
        else
            testT1map(:) = T1;
            testT2map(:) = T2;
        end
        
        % add some small variation in the starting T1 and T2 values
        testT1map = testT1map.*imnoise(testT1map,'gaussian',0, (0.005)^2); %Fourth argument is the variance
        testT2map = testT2map.*imnoise(testT2map,'gaussian',0, (0.005)^2);
        

        
        
end
end