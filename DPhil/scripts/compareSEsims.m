% script to test three simulation functions

% parameters
T1 = 1000;
 T2 = 100;
 nPts = 96;

 %make offset list
offsets = generate_offset_list([130 32 75 165],[5000 82 105 195],nPts,[1 2 3 4],'SSEPI','nosave');
 
[sigEPG, P] = sim_SE_EPG(complex(deg2rad(offsets(:,3)),0),complex(deg2rad(offsets(:,4)),0), offsets(:,1)/1000, offsets(:,2)/1000, T1/1000, T2/1000);
sigB = sim_SE_bernstein(T1, T2, 0, 0, offsets,1,0);
sigH = sim_SE_Hargreaves(T1, T2, 0, 0, offsets,1, 0, 'noplot');

%plot results
figure
 plot(abs(sigEPG/norm(sigEPG)),'g-o')
 hold on
 plot(abs(sigB/norm(sigB)),'b-*')
 plot(abs(sigH/norm(sigH)),'r-+')
 legend('EPG','Bernstein','Hargeaves')
 xlabel 'Time Point'
 ylabel 'Normalised Absolute Signal'