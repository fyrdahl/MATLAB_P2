% IRbSSFP simulation figure

for n=1:4
    for m = 1:size(IRbSSFPT1err,3)
    
    T1stds(n,m) = std(squeeze(IRbSSFPT1err(n,1,m,1,:)))
    
    end
    
end

figure, 
x = [50,100:100:1000];
range = 1:size(IRbSSFPT1err,3)
errorbar(x,mean(squeeze(IRbSSFPT1err(1,1,range,1,:)),2),squeeze(T1stds(1,range)))
hold on
errorbar(x,mean(squeeze(IRbSSFPT1err(2,1,range,1,:)),2),squeeze(T1stds(2,range)))
errorbar(x,mean(squeeze(IRbSSFPT1err(3,1,range,1,:)),2),squeeze(T1stds(3,range)))
errorbar(x,mean(squeeze(IRbSSFPT1err(4,1,range,1,:)),2),squeeze(T1stds(4,range)))


for n=1:4
    for m = 1:12
    
    T2stds(n,m) = std(squeeze(IRbSSFPT2err(n,1,m,1,:)));
    
    end
    
end

figure, 
x = [50,100:100:1000];
range = 2:12
errorbar(x,mean(squeeze(IRbSSFPT2err(1,1,range,1,:)),2),squeeze(T2stds(1,range)))
hold on
errorbar(x,mean(squeeze(IRbSSFPT2err(2,1,range,1,:)),2),squeeze(T2stds(2,range)))
errorbar(x,mean(squeeze(IRbSSFPT2err(3,1,range,1,:)),2),squeeze(T2stds(3,range)))
errorbar(x,mean(squeeze(IRbSSFPT2err(4,1,range,1,:)),2),squeeze(T2stds(4,range)))