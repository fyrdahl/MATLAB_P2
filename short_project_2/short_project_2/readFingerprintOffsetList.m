%% read fingerprint offset timings list

fileID = fopen('Fingerprint_List.txt','r');
fileID = fopen('/Users/jallen/Documents/short_project_2/20150815/Fingerprint_List.txt','r');
(fscanf(fileID,'%s',2));
nList1 = str2num((fscanf(fileID,'%s',1)));
for i = 1:nList1
    fingerprintLists(i,:,1) = fscanf(fileID,'%f',4);
end
(fscanf(fileID,'%s',2));
nList2 = str2num((fscanf(fileID,'%s',1)));
for i = 1:nList2
    fingerprintLists(i,:,2) = fscanf(fileID,'%f',4);
end
(fscanf(fileID,'%s',2));
nList3 = str2num((fscanf(fileID,'%s',1)));
for i = 1:nList3
    fingerprintLists(i,:,3) = fscanf(fileID,'%f',4);
end
(fscanf(fileID,'%s',2));
nList4 = str2num((fscanf(fileID,'%s',1)));
for i = 1:nList4
    fingerprintLists(i,:,4) = fscanf(fileID,'%f',4);
end
(fscanf(fileID,'%s',2));
nList5 = str2num((fscanf(fileID,'%s',1)));
for i = 1:nList5
    fingerprintLists(i,:,5) = fscanf(fileID,'%f',4);
end
(fscanf(fileID,'%s',2));
nList6 = str2num((fscanf(fileID,'%s',1)));
for i = 1:nList6
    fingerprintLists(i,:,6) = fscanf(fileID,'%f',4);
end
(fscanf(fileID,'%s',2));
nList7 = str2num((fscanf(fileID,'%s',1)));
for i = 1:nList7
    fingerprintLists(i,:,7) = fscanf(fileID,'%f',4);
end
(fscanf(fileID,'%s',2));
nList8 = str2num((fscanf(fileID,'%s',1)));
for i = 1:nList8
    fingerprintLists(i,:,8) = fscanf(fileID,'%f',4);
end
(fscanf(fileID,'%s',2));
nList9 = str2num((fscanf(fileID,'%s',1)));
for i = 1:nList9
    fingerprintLists(i,:,9) = fscanf(fileID,'%f',4);
end
fclose(fileID);

