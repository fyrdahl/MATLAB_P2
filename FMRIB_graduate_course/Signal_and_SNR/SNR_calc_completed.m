function [snr] = SNR_calc_completed(image,sig_mask,bg_mask)

% this is the region we are calculating the SNR in
sig_index=find(sig_mask==1);
% this is our background region where we estimate the noise
bg_index=find(bg_mask==1);
% calculate the mean signal inside the signal mask
sig_mean=mean(image(sig_index));
% TASK 1: calculate the standard deviation inside the background mask
bg_std=std(image(bg_index));
% TASK 2: calculate the SNR (don't forget to multiply by 0.65!)
snr=0.65*sig_mean/bg_std;
% TASK 3: print out the results
disp(['mean signal: ' num2str(sig_mean)]); disp(['noise std: ' num2str(bg_std)]); disp(['SNR: ' num2str(snr)]);

  
  