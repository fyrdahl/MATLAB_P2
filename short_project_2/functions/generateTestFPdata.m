function testFPtcMap = generateTestFPdata(nTimePts,fingerprintOffsetList,freqOffset,nSlices,workingdir)
%%
% nTimePts: Defines the number points in the generated timecourse
%
%%
sIx = 20;
sIy = 20;

readFingerprintOffsetList;
%%
testFPtcMap = zeros(sIx*sIy,nTimePts);
M = zeros(1,nTimePts);
M0=1;
%%
testT1 = zeros(sIx,sIy);
testT1(:) = 220;
testT1(10:15,10:15) = 260;
testT1(5:10,5:10) = 280;
testT2 = testT1*0.8;
testT1 = reshape(testT1,[1,sIx*sIy]);
testT2 = reshape(testT2,[1,sIx*sIy]);
%%
for i = 1:size(testT1,2)
 
    [~, Mxy] = SimSE_Bernstein(testT1(i), testT2(i), fingerprintOffsetList,freqOffset, nSlices, nTimePts);
    testFPtcMap(i,:) = Mxy;
end
    testFPtcMap = reshape(testFPtcMap,[sIx,sIy,nTimePts]);

end
