offsetList = generate_offset_list([900, 100, 70, 180],[1100 100 90 180],24,1);
testFPtcMap = generate_test_FP_data(24,offsetList,0,1);
[M, Mxy,flipAngles,imageTimes, t0s] = sim_SE_bernstein(T1, T2, offsetList,0, 1, 24,'plot');

