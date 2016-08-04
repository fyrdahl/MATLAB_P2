<<<<<<< HEAD
%% Script to test the simulation of a chosen offset list

T1 = 282.3;
T2 = 214.8;
nTimePts = 24;
% offsetList = generate_offset_list(lower,upper,nPts,seed,seq)
offsetList = generate_offset_list([900, 100, 90, 180],[900 100 90 180],nTimePts,1,'SE_EPI');
% [testFPtcMap, testT1map, testT2map] = generate_test_FP_data(nTimePts,offsetList,freqOffset,nSlices)
testFPtcMap = generate_test_FP_data(24,offsetList,0,1);
%[M, Mxy,flipAngles,imageTimes, t0s] = sim_SE_bernstein(T1, T2, offsetList,0, 1, 24,'plot');
%[M, Mxy,imageTimes,flipAngles, t0s] = sim_SE(T1, T2, fingerprintOffsetList,nRepeats,df, nSlices, plotFlag)
[M, Mxy,imageTimes,flipAngles, t0s] = sim_SE(T1, T2, offsetList,(nTimePts/24),0, 1, 'plot');
=======
offsetList = generate_offset_list([900, 100, 70, 180],[1100 100 90 180],24,1);
testFPtcMap = generate_test_FP_data(24,offsetList,0,1);
[M, Mxy,flipAngles,imageTimes, t0s] = sim_SE_bernstein(T1, T2, offsetList,0, 1, 24,'plot');

>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
