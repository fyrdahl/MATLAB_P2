function offsetList = generate_offset_list(lower,upper,nPts,seed)
% GENERATE_OFFSET_LIST  Offset list for use with spin echo fingerprinting
%
% Function to generate a list of parameters to be used when running spin
% echo sequence with variable TR, TE, excitation flip angle (FA1) and
% refocusing flip angle (FA2).
%
% [OFFSETLIST] = GENERATE_OFFSET_LIST(LOWER, UPPER, NPTS, SEED)
%
%   LOWER is 4-element vector of the lower limits for a uniform distribution
%   of each parameter (TR, TE, FA1, FA2).
%
%   UPPER is the upper limits.
%
%   NPTS is the number of timepoints wanted from the sequence.
%
%   SEED determines which random numbers are generated.

% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
% 2015

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

end

