function offsetList = generate_offset_list(lower,upper,nPts,seed,seq,saveFlag,varargin)
% GENERATE_OFFSET_LIST  Offset list for use with spin echo fingerprinting
%
% Function to generate a list of parameters to be used when running spin
% echo sequence with variable TR, TE, excitation flip angle (FA1) and
% refocusing flip angle (FA2).
%
% [OFFSETLIST] = GENERATE_OFFSET_LIST(LOWER, UPPER, NPTS, SEED)
%
%   LOWER is 4-element vector of the lower limits for a uniform distribution
%   of each parameter (TR, TE, FA1, FA2) in micro-seconds.
%
%   UPPER is the upper limits in micro-seconds.
%
%   NPTS is the number of timepoints wanted from the sequence.
%
%   SEED is a vector, determines which random numbers are generated.

% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
% 2015


switch seq
    case 'SSEPI'
        % create a list
        offsetList = zeros(nPts,numel(lower));
        rng(seed(1),'v5uniform')
        offsetList(:,1) = round((upper(1)-lower(1)).*rand(nPts,1) + lower(1)); %TRs in ms
        rng(seed(2),'v5uniform')
        offsetList(:,2) = round((upper(2)-lower(2)).*rand(nPts,1) + lower(2));%TEs in ms
        rng(seed(3),'v5uniform')
        offsetList(:,3) = round((upper(3)-lower(3)).*rand(nPts,1) + lower(3));%FA1s
        rng(seed(4),'v5uniform')
        offsetList(:,4) = round((upper(4)-lower(4)).*rand(nPts,1) + lower(4));%FA2s
        
        if strcmp(saveFlag,'save')
            %add the list to end of 'Fingerprint_List.txt'
            fID = fopen(['~/Documents/MATLAB/DPhil/offsetLists/Fingerprint_List.txt'],'a');
            fprintf(fID,['LIST_OFFSET ',num2str(varargin{1}),'\r\n']);
            fprintf(fID,[num2str(nPts),'\r\n']);
            for i = 1:size(offsetList,1)
                for k = 1:3
                    fprintf(fID,[num2str(offsetList(i,k)),' ']);
                end
                fprintf(fID,num2str(offsetList(i,4)));
                fprintf(fID,'\r\n');
            end
            disp('offset list appended to Fingerprint_List.txt')
            
        end
        
    case 'SSFP'
        % create a list
        offsetList = zeros(nPts,numel(lower));
        rng(seed(1),'v5uniform')
        offsetList(:,1) = (upper(1)-lower(1)).*rand(nPts,1) + lower(1); %TRs in ms
        rng(seed(2),'v5uniform')
        offsetList(:,2) = (upper(2)-lower(2)).*rand(nPts,1) + lower(2);%FAs
        
        if strcmp(saveFlag,'save')
            %add the list to end of 'Fingerprint_List.txt'
            fIDssfp = fopen(['~/Documents/MATLAB/DPhil/offsetLists/NEW_SSFP_Fingerprint_List.txt'],'a');
            fprintf(fIDssfp,['LIST_OFFSET ',num2str(varargin{1}),'\r\n']);
            fprintf(fIDssfp,[num2str(nPts),'\r\n']);
            for i = 1:size(offsetList,1)
                    fprintf(fIDssfp,[num2str(offsetList(i,1)),' ']);
                    fprintf(fIDssfp,num2str(offsetList(i,2)));
                fprintf(fIDssfp,'\r\n');
            end
            disp('offset list appended to Fingerprint_List.txt')
        end
end



end

