function [f, absMsig1, absMsig2] = dot_product_sim_SE(params,relaxParams)

nPts = 48;
nRepeats =  nPts/24;
offsetList = zeros(48,4);
offsetList(:,1) = repmat(params(1),nPts,1);
offsetList(:,2) = repmat(params(2),nPts,1);
offsetList(:,3) = repmat(params(3),nPts,1);
offsetList(:,4) = repmat(params(4),nPts,1);
df = 0;

[absMsig1] = sim_SE3(relaxParams(1), relaxParams(2), offsetList,nRepeats,df, 'noPlot');
[absMsig2] = sim_SE3(relaxParams(3), relaxParams(4), offsetList,nRepeats,df, 'noPlot');
absMsig1
absMsig2
f = sum(abs(absMsig1-absMsig2));
end

