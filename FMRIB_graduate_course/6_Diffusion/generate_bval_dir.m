function [bvalall,diffdirall] = generate_bval_dir(handles)

% find a bug. the siemens diffusion file in not normalised


dirbegin = ['[directions=',num2str(handles.nodiffdir),']'];
[fid,errstr]=fopen('DiffusionVectors.txt');
if ( fid == -1 ),
    error(errstr);
end;

i_loop_counter = 0;
tag1 = 0;
diffdirall = zeros(3,handles.nodiffdir);
bvalall = ones(1,handles.nodiffdir)*handles.bvals;
ttt = 0;
ispe = 0;
while 1
    i_loop_counter = i_loop_counter + 1;
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    if strcmp(tline, dirbegin)
    tag1 = 1;
    ispe = i_loop_counter;

    end
    
    if (i_loop_counter>= (ispe+3))&&(i_loop_counter<=(ispe + handles.nodiffdir + 2 )&&(tag1 == 1))
      ttt = ttt  + 1 ;
       ppp1 = findstr(tline,'(');
       ppp2 = findstr(tline,')');
       str = tline((ppp1+1) : (ppp2-1));
       dirtmp = textscan(str,'%f','Delimiter',','); 
       dirtmp = cell2mat(dirtmp);
       diffdirall(:,ttt) = dirtmp;
    end
end

% normalization
for i = 1 : handles.nodiffdir

    normalvalue =  sqrt(diffdirall(1,i)^2 + diffdirall(2,i)^2 + diffdirall(3,i)^2);
    if normalvalue ~=0
    diffdirall(:,i) = diffdirall(:,i)/normalvalue;
    end

end

fclose(fid);

fid = fopen('bvecs','w');

for ii = 1 : 3
    for jj = 1 : handles.nodiffdir
        if diffdirall(ii,jj) >= 0
            fprintf(fid,'   %6.6f',diffdirall(ii,jj));
        else
            fprintf(fid,'  %6.6f',diffdirall(ii,jj));
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);


for ii = 1 : handles.nodiffdir
   
    if sum(abs(diffdirall(:,ii)),1) == 0
        bvalall(ii) = 0;
    end
    
end


fid = fopen('bvals','w');

    for jj = 1 : handles.nodiffdir
       
            fprintf(fid,'  %d',bvalall(jj));


    end

fclose(fid);









