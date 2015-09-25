%
switch phantomName
    case 'Jack'
        
        %TE
        compartmentCenters(1,1:2,1) = [36 23];
        compartmentCenters(2,1:2,1) = [28 21];
        compartmentCenters(3,1:2,1) = [22 28];
        compartmentCenters(4,1:2,1) = [25 36];
        compartmentCenters(5,1:2,1) = [34 38];
        compartmentCenters(6,1:2,1) = [39 31];
        
        %TI
        compartmentCenters(1,1:2,2) = [compartmentCenters(1,1:2,1) + 3];
        compartmentCenters(2,1:2,2) = [compartmentCenters(2,1:2,1) + 3];
        compartmentCenters(3,1:2,2) = [compartmentCenters(3,1:2,1) + 3];
        compartmentCenters(4,1:2,2) = [compartmentCenters(4,1:2,1) + 3];
        compartmentCenters(5,1:2,2) = [compartmentCenters(5,1:2,1) + 3];
        compartmentCenters(6,1:2,2) = [compartmentCenters(6,1:2,1) + 3];
        
        % fingerprinting
        compartmentCenters(1,1:2,3) = [37 23];
        compartmentCenters(2,1:2,3) = [28 24];
        compartmentCenters(3,1:2,3) = [25 32];
        compartmentCenters(4,1:2,3) = [30 38];
        compartmentCenters(5,1:2,3) = [38 37];
        compartmentCenters(6,1:2,3) = [42 30];
        
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
      
        compartmentCenters(1,1:2,3) = [21 21];
        compartmentCenters(2,1:2,3) = [22 22];
        compartmentCenters(3,1:2,3) = [24 24];
        compartmentCenters(4,1:2,3) = [27 27];
        compartmentCenters(5,1:2,3) = [32 32];
        compartmentCenters(6,1:2,3) = [39 42];
      
        
    
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


