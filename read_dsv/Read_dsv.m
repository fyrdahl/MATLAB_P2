function dsv = Read_dsv (dsv_filename )
% Read Siemens dsv
%
%
% Read spectroscopy data from Siemens machine
%
% Read a .dsv file
% put the tags into the dsv structure
% store the data into dsv.output_data
%
% Then convert on the output into the right form (uncompressing
% it) (using dsv2timecourse) and output into dsv.timecourse
%
%
% If filename is not supplied, then ask for one
%
% MDR April 2008
%
if (nargin < 1)
    [input_filename, input_pathname] = uigetfile ('*.dsv','Pick one of the input files');
    dsv_filename = [input_pathname '/' input_filename];
end

dsv.filename = dsv_filename;



%dsv_filename = 'c:/Temp/DspData_M0Y.dsv';

fid = fopen(dsv_filename);

comment_text = ';';
head_end_text   = '[VALUES]';

tline = fgets(fid);

while (isempty(strfind(tline , head_end_text)))
    
    tline = fgets(fid);
    
    if ( strcmp(tline(1) , comment_text) == 0 & size(deblank(tline)) > 0 )
        
        
        % Store this data in the appropriate format
        
        occurence_of_colon = findstr('=',tline);
        variable = tline(1:occurence_of_colon-1) ;
        value    = tline(occurence_of_colon+1 : length(tline)) ;
        
        if (~isempty(variable))
        switch variable
            case { 'FORMAT' ,  'TITLE' , 'VERTUNITNAME' , 'HORIUNITNAME'  }
                eval(['dsv.' , variable , ' = value ;']);
                
                
            case {  'SAMPLES'  }
                %Integers
                eval(['dsv.' , variable , ' = str2num(value) ;']);
            case { 'VERTFACTOR' , 'MINLIMIT' , 'MAXLIMIT' , 'HORIDELTA' }
                %Floats 
                eval(['dsv.' , variable , ' = str2num(value) ;']);
                
            case { ''  }
                
            otherwise
                % We don't know what this variable is.  Report this just to keep things clear
                disp(['Unrecognised variable ' , variable ]);
        end
    end
        %else
        % Don't bother storing this bit of the output
        %end
        
    end
end


% Now read the data line by line until we get to the bit at the end of the
% file, where it says that there is no more information!

%

output_data_counter = 1;
tline = fgets(fid);
while (  tline  ~= -1 );
    
    if ( strcmp(tline(1) , comment_text) == 0 & size(deblank(tline)) > 0 )
            
        % Store this data in the appropriate format
        value    = str2num(tline);
        dsv.output_data(output_data_counter) = value;
        output_data_counter = output_data_counter + 1;
        
    end
    tline = fgets(fid);
end


dsv.timecourse = dsv2timecourse (dsv);

fclose(fid);