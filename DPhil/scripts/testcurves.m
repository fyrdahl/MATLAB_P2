%% plot test magnestisation evolution curves
Mzeq = 1 ;
Mz0 = -Mzeq ;
B = 1 ;
x = [0:10:1500] ;
T1list = [250; 193; 143; 126; 52] ;
T2list = [98; 71; 53; 43; 21] ;



x = TI
for i = T1list'
   
    s = Mzeq - (Mzeq - Mz0)*B*exp(-x*(1/i));
      
%  noise = (noiseStd*randn(1,numel(s))) + noiseMean;
%  noise = noise/signalMean;
%      y(i,:) = s + noise;
%      
% [fittedCurve, goodness, output] = fit(x',y(i,:)',T1model,'Upper',[2000 5 10 5000],'Lower',[0 -5 0 0],'StartPoint',[1000 100 1 100])
%    sse(i) = goodness.sse ;
end


figure; hold on
plot(s)
%%
plot(x, y,'o')
plot(x,s,'g-')
plot(fittedCurve,'r-')


%%
figure;
for T2 = T2list'
    hold on
plot(x, Mzeq*exp(-x*(1/T2)))
end
legend (num2str(T2list))