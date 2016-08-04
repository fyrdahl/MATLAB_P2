<<<<<<< HEAD
function [syntheticFPdata] = generate_test_FP_data(testT1map, testT2map, minTR, minTE, nRuns,offsets,df,varargin)
=======
function testFPtcMap = generate_test_FP_data(nTimePts,offsetList,freqOffset,nSlices)
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
%% Generate Fingerprint Time Course using synthesised maps of T1 and T2
% nTimePts: Defines the number points in the generated timecourse
%
%
<<<<<<< HEAD
tic
disp(['Test T1 map dimensions: ',num2str(size(testT1map))])
disp(['Test T2 map dimensions: ',num2str(size(testT2map))])
%% Simulate fingerprint signal using test T1 and T2 maps
xdim = size(testT1map,2);
ydim = size(testT1map,1);
nTimeCoursePts = nRuns*size(offsets,1);
syntheticFPdata = zeros(nTimeCoursePts, xdim*ydim);

for nVoxel = 1:xdim*ydim
    
    if (testT1map(nVoxel)>0) && (testT2map(nVoxel)>0)
        
        
        % function [Msig, M] = sim_SE_bernstein(T1, T2, minTR, minTE, offsets,nRepeats,df)
        if nargin > 8
            B1 = varargin{2};            
            tmpoffsets(:,1) = offsets(:,1);
            tmpoffsets(:,2) = offsets(:,2);
            tmpoffsets(:,3) = offsets(:,3)*B1(nVoxel);
            tmpoffsets(:,4) = offsets(:,4)*B1(nVoxel);
            
        else
            tmpoffsets = offsets;
        end
       
        
        
        syntheticFPdata(:,nVoxel) = sim_SE_bernstein(testT1map(nVoxel), testT2map(nVoxel), minTR, minTE, tmpoffsets, nRuns,df);
        
        
        if nargin > 7
            % add white gaussian noise to each synthetic fp timecourse
            syntheticFPdata(:,nVoxel) = awgn(syntheticFPdata(:,nVoxel),varargin{1});
        end
        
    end
end

toc
=======
%% Make Test T1 and T2 maps
sIx = 20;
sIy = 20;
testFPtcMap = zeros(sIx*sIy,nTimePts);
testT1 = zeros(sIx,sIy);
testT1(:) = 150;
testT1(10:15,10:15) = 250;
testT1(5:10,5:10) = 300;
testT2 = testT1*0.9;
testT1 = reshape(testT1,[1,sIx*sIy]);
testT2 = reshape(testT2,[1,sIx*sIy]);

%% Simulate fingerprint signal using test T1 and T2 maps
for i = 1:size(testT1,2)
    [~, Mxy] = sim_SE_bernstein(testT1(i), testT2(i), offsetList,freqOffset, nSlices, nTimePts,'no plot');
    testFPtcMap(i,:) = Mxy;
end
%Transfer back into 3D matrix form
testFPtcMap = reshape(testFPtcMap,[sIx,sIy,nTimePts]);
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
end
