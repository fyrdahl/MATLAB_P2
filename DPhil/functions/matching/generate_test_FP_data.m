function [syntheticFPdata] = generate_test_FP_data(testT1map, testT2map, minTR, minTE, nRuns,offsets,df,varargin)
%% Generate Fingerprint Time Course using synthesised maps of T1 and T2
% nTimePts: Defines the number points in the generated timecourse
%
%
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
end
