%% generate optimal offset list

% range of properties to use to find the optimal list
T1 = [ 300 ; 250 ];
T2 = [ 250; 200 ];

% create an offsetList from which to start
nPts = 5;
df = 0;
nRepeats = 1;
nSlices = 1;
lower = [400 0 75 165]; % TR, TE, FA1, FA2
upper = [3000 85 105 195];
TEmin = 32; %ms
TRmin = 130;

TRadd = [ 1:10]*20;
TEadd = [ 1]*2000;
FA1add = [ 1]*200;
FA2add = [ 1]*200;

%key: offsetList = generate_offset_list(lower,upper,nPts,seed,seq,offsetNumber,saveFlag)
offsetList = generate_offset_list(lower,upper,nPts,[1 2 3 4],'SE-EPI','dontSave');

for indTR = 1:numel(TRadd)
    for indTE = 1:numel(TEadd)
        for indFA1 = 1:numel(FA1add)
            for indFA2 = 1:1:numel(FA2add)
                for n = 1:numel(T1)
                    updatedOffsetList(:,1) = offsetList(:,1)+TRadd(indTR);
                    updatedOffsetList(:,2) = offsetList(:,2)+TEadd(indTE);
                    updatedOffsetList(:,3) = offsetList(:,3)+FA1add(indFA1);
                    updatedOffsetList(:,4) = offsetList(:,4)+FA2add(indFA2);
                    simTC(n,:) = sim_SE3(T1(n), T2(n), TRmin,TEmin,updatedOffsetList,nRepeats,1, df, 'dontPlot');
                end
                NDP(indTR,indTE,indFA1,indFA2) = dot(simTC(1,:),simTC(2,:))/(norm(simTC(1,:))*norm(simTC(2,:)))
                updatedOffsetList
                pause
            end
        end
    end
    
end
min(NDP(:))
NDP(find(NDP == (min(NDP(:)))))
NDP

%%
optTE = t(var(signal) == max(var(signal)));

plot(Mz)
offsetList = generate_offset_list([max(t(Mz <= 0.05*1)) 12 90 180],[min(t(Mz >= 0.95*1)) 645 90 180],24,1,'SE-EPI');
plot(round(nT2),signal(round(nT2)),'*')


freqOffset = 0;
nSlices = 2;
nRepeats = 2;
nTimeCoursePts = nRepeats*size(offsetList,1);
TRs = 10:100;

for nTR = 1:10
    
    for nT2 = 1:2
        T2s = [100,260];
        [absMsig(nT2)] = sim_SE3(T1, T2(nT2), offsetList,nRepeats,df, plotFlag);
    end
    sigDiff(nTR) = sum(absMsig(1,:) - absMsig(2,:));
end












