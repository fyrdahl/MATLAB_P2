% 
% compartmentCenters(1,1:2) = [37 24]
% compartmentCenters(2,1:2) = [28 25]
% compartmentCenters(3,1:2) = [25 32]
% compartmentCenters(4,1:2) = [30 39]
% compartmentCenters(5,1:2) = [39 38]
% compartmentCenters(6,1:2) = [42 30]

% compartmentCenters(1,1:2) = [28 22]
compartmentCenters(1,1:2) = [21 21]
% compartmentCenters(2,1:2) = [28 25]
compartmentCenters(2,1:2) = [22 22]
% compartmentCenters(3,1:2) = [25 32]
compartmentCenters(3,1:2) = [23 23]
% compartmentCenters(4,1:2) = [30 39]
compartmentCenters(4,1:2) = [41 21]
compartmentCenters(5,1:2) = [39 38]
compartmentCenters(6,1:2) = [42 30]

TE = TE
figure
for i = 1:6
    hold on
plot(TE(2:end),squeeze(TEimages(compartmentCenters(i,1),compartmentCenters(i,2),TE(2:end))),'-*')

end

 TI = TI
figure
for i = 1:6
    hold on
    compartmentTIimages(i,:) = TIimages(compartmentCenters(i,1),compartmentCenters(i,2),TI(2:11))
plot(TI(2:11),compartmentTIimages(i,:),'-*')

end


