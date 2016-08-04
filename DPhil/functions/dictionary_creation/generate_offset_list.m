<<<<<<< HEAD
function offsetList = generate_offset_list(lower,upper,nPts,seed,seq,saveFlag,varargin)
=======
function offsetList = generate_offset_list(lower,upper,nPts,seed)
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
% GENERATE_OFFSET_LIST  Offset list for use with spin echo fingerprinting
%
% Function to generate a list of parameters to be used when running spin
% echo sequence with variable TR, TE, excitation flip angle (FA1) and
% refocusing flip angle (FA2).
%
% [OFFSETLIST] = GENERATE_OFFSET_LIST(LOWER, UPPER, NPTS, SEED)
%
%   LOWER is 4-element vector of the lower limits for a uniform distribution
<<<<<<< HEAD
%   of each parameter (TR, TE, FA1, FA2) in micro-seconds.
%
%   UPPER is the upper limits in micro-seconds.
%
%   NPTS is the number of timepoints wanted from the sequence.
%
%   SEED is a vector, determines which random numbers are generated.
=======
%   of each parameter (TR, TE, FA1, FA2).
%
%   UPPER is the upper limits.
%
%   NPTS is the number of timepoints wanted from the sequence.
%
%   SEED determines which random numbers are generated.
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1

% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
% 2015

<<<<<<< HEAD

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



=======
% create a list
offsetList = zeros(nPts,numel(lower));
rng(seed,'twister')
offsetList(:,1) = round(unifrnd(lower(1),upper(1),nPts,1)); %TRs
rng(seed,'twister')
offsetList(:,2) = round(unifrnd(lower(2),upper(2),nPts,1));%TEs
rng(seed,'twister')
offsetList(:,3) = round(unifrnd(lower(3),upper(3),nPts,1));%FA1s
rng(seed,'twister')
offsetList(:,4) = round(unifrnd(lower(4),upper(4),nPts,1));%FA2s

% save the newly generated list, with a unique file name (using the date
% and time).
date = datestr(now,30);
fid = fopen(['/Users/jallen/Documents/MATLAB/DPhil/offsetLists/',date,'Fingerprint_List'],'w');
fprintf(fid,'LIST_OFFSET 1 \n');
for i = 1:size(offsetList,1)
    for k = 1:4
fprintf(fid,[num2str(offsetList(i,k)),' ']);
    end
fprintf(fid,'\n');
end

>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
end

