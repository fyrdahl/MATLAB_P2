% IRbSSFP simulation figure

for n=1:4
    for m = 1:size(SSEPIT1err,3)
    
    SSEPIT1errstds(n,m) = std(squeeze(SSEPIT1err(n,1,m,1,:)));
    
    end
    
end

figure, 

range = 1:size(SSEPIT1err,3)
errorbar(mean(squeeze(SSEPIT1err(1,1,range,1,:)),2),squeeze(SSEPIT1errstds(1,range)))
hold on
errorbar(mean(squeeze(SSEPIT1err(2,1,range,1,:)),2),squeeze(SSEPIT1errstds(2,range)))
errorbar(mean(squeeze(SSEPIT1err(3,1,range,1,:)),2),squeeze(SSEPIT1errstds(3,range)))
errorbar(mean(squeeze(SSEPIT1err(4,1,range,1,:)),2),squeeze(SSEPIT1errstds(4,range)))


for n=1:4
    for m = 1:size(SSEPIT2err,3)
    
    SSEPIT2errstds(n,m) = std(squeeze(SSEPIT2err(n,1,m,1,:)));
    
    end
    
end

figure, 
range = 1:size(SSEPIT1err,3)
errorbar(mean(squeeze(SSEPIT2err(1,1,range,1,:)),2),squeeze(SSEPIT2err(1,range)))
hold on
errorbar(mean(squeeze(SSEPIT2err(2,1,range,1,:)),2),squeeze(SSEPIT2err(2,range)))
errorbar(mean(squeeze(SSEPIT2err(3,1,range,1,:)),2),squeeze(SSEPIT2err(3,range)))
errorbar(mean(squeeze(SSEPIT2err(4,1,range,1,:)),2),squeeze(SSEPIT2err(4,range)))