function testFPtcMap = generateTestFPdata(nTimePts,offsetListNum,freqOffset,nSlices,workingdir)
%% Generate Fingerprint Time Course using synthesised maps of T1 and T2
% nTimePts: Defines the number points in the generated timecourse
%
%
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
offsetList = readFPOffsetList(workingdir,offsetListNum);
for i = 1:size(testT1,2)
    [~, Mxy] = SimSE_Bernstein(testT1(i), testT2(i), offsetList,freqOffset, nSlices, nTimePts);
    testFPtcMap(i,:) = Mxy;
end
%Transfer back into 3D matrix form
testFPtcMap = reshape(testFPtcMap,[sIx,sIy,nTimePts]);
end
