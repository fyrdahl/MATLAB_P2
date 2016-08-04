function SNR_calc_with_ROIs(image,max_intensity)

if nargin < 2; max_intensity = max(image(:)); end

disp('First draw the signal ROI')
sig_mask = make_roi(image,'First draw the signal ROI');

disp('Now draw the background ROI')
bg_mask = make_roi(image,'Now draw the background ROI',max_intensity);

disp('Calculating SNR')
SNR_calc_completed(image,sig_mask,bg_mask);