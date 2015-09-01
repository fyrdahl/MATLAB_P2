%
switch phantomName
    case 'Jack'
        
        %TE
        compartmentCenters(1,1:2,1) = [23 36];
        compartmentCenters(2,1:2,1) = [21 28];
        compartmentCenters(3,1:2,1) = [28 22];
        compartmentCenters(4,1:2,1) = [36 25];
        compartmentCenters(5,1:2,1) = [38 34];
        compartmentCenters(6,1:2,1) = [31 39];
        
        %TI
        compartmentCenters(1,1:2,2) = [compartmentCenters(1,1:2,1) + 3];
        compartmentCenters(2,1:2,2) = [compartmentCenters(2,1:2,1) + 3];
        compartmentCenters(3,1:2,2) = [compartmentCenters(3,1:2,1) + 3];
        compartmentCenters(4,1:2,2) = [compartmentCenters(4,1:2,1) + 3];
        compartmentCenters(5,1:2,2) = [compartmentCenters(5,1:2,1) + 3];
        compartmentCenters(6,1:2,2) = [compartmentCenters(6,1:2,1) + 3];
        
        % fingerprinting
        compartmentCenters(1,1:2,3) = [23 37];
        compartmentCenters(2,1:2,3) = [24 28];
        compartmentCenters(3,1:2,3) = [32 25];
        compartmentCenters(4,1:2,3) = [38 30];
        compartmentCenters(5,1:2,3) = [37 38];
        compartmentCenters(6,1:2,3) = [30 42];
        
    case 'sphereD170'
       
        compartmentCenters(1,1:2,1) = [21 21];
        compartmentCenters(2,1:2,1) = [22 22];
        compartmentCenters(3,1:2,1) = [24 24];
        compartmentCenters(4,1:2,1) = [27 27];
        compartmentCenters(5,1:2,1) = [32 32];
        compartmentCenters(6,1:2,1) = [39 42];
        
        compartmentCenters(1,1:2,2) = [21 21];
        compartmentCenters(2,1:2,2) = [22 22];
        compartmentCenters(3,1:2,2) = [24 24];
        compartmentCenters(4,1:2,2) = [27 27];
        compartmentCenters(5,1:2,2) = [32 32];
        compartmentCenters(6,1:2,2) = [39 42];
      
        
    
end


figure
for i = 1:6
    hold on
    plot(TE(2:end),squeeze(TEimages(compartmentCenters(i,1,1),compartmentCenters(i,2,1),TE(2:end))),'-*')
    
end


figure
for i = 1:6
    hold on
    plot(TI(2:end), squeeze(TIimages(compartmentCenters(i,1,2),compartmentCenters(i,2,2),TI(2:end))),'-*')
  
end


