function matchout=templatematch_sim(dict,sig)
%
% This function calculates a template match for use in MR Fingerprinting
% applications. Here this is calculated using a simple dot-product. The
% output is the index of the most likely dictionary entry for each pixel.
%
% INPUTS: dict (# entries, timepoints)
% sig (timepoints, #pixels)
%
% OUTPUTS: matchout (1, #pixels)
%
% Please note that this code requires the generation of an intermediate
% matrix of size (#entries, #pixels), which can be quite large. If this is
% beyond the memory capability of your computer, simply parse the data into
% smaller sets of pixels and calculate them serially.
%
% Nov. 15, 2012 Dan Ma, Mark Griswold
%
if nargin<2
disp('Error - 2 inputs are required!');
return
end
% First normalize the dictionary and signal
for i=1:size(dict(:,:),1)
dict(i,:)=dict(i,:)./norm(dict(i,:)); %JA: fixed bug in original code (was not normalised in correct direction)
end
for i=1:size(sig(:,:),2)
sig(:,i)=sig(:,i)./norm(sig(:,i));
end
% Calculate inner product
innerproduct=dict(:,:)*sig(:,:);
% Take the maximum value and return the index
[maxval,matchout]=max(abs(innerproduct));

end