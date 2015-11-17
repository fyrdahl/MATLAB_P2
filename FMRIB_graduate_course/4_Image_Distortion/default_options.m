function options=default_options(options)
  
  if (nargin<1) 
      options=struct(); 
  else 
      orig_fields = fieldnames(orderfields(options));
      recognised_fields = ({'GRAPPA'
          'AddNoise';
          'Bandwidth';
          'FatSat';
          'FieldStrength';
          'BlipUp';
          'Matrix';
          'Dropout';
          'SliceThickness';
          'AddSinus';
          'TE';
          'PhaseEncodeLR';});
      if any(~ismember((orig_fields), recognised_fields))
          indx = ~ismember((orig_fields), recognised_fields);
          disp(['*********************************************************'])
          error(['*** I do not recognise options.' cell2mat(orig_fields(find(indx,1,'first'))) '. You can remove this field by using:     options=rmfield(options,''' cell2mat(orig_fields(find(indx,1,'first'))) ''')'])    
 
      end         
  end;
          

  % first, catch capitalization mistakes
  if isfield(options,'grappa')  
    options.GRAPPA=options.grappa;  end;
  if isfield(options,'bandwidth')
    options.Bandwidth=options.bandwidth; end;
  if isfield(options,'fieldstrength')  
    options.FieldStrength=options.fieldstrength; end;
  if isfield(options,'matrix')         
    options.Matrix=options.matrix; end;
  if isfield(options,'slicethickness')
    options.SliceThickness= options.slicethickness; end;
  if isfield(options,'te')             
    options.TE=options.te; end;

  if isfield(options,'addnoise')       
    options.AddNoise=options.addnoise; end;
  if isfield(options,'fatsat')         
    options.FatSat=options.fatsat; end;
  if isfield(options,'blipup')         
    options.BlipUp=options.blipup; end;
  if isfield(options,'dropout')        
    options.Dropout=options.dropout; end;
  if isfield(options,'addsinus')       
    options.AddSinus=options.addsinus; end;
  if isfield(options,'phaseencodelr')
    options.PhaseEncodeLR=options.phaseencodelr; end;
  if isfield(options,'phaseencodeLR')
    options.PhaseEncodeLR=options.phaseencodeLR; end;

    
  % set uninitialized options to default
  if ~isfield(options,'GRAPPA')         
    options.GRAPPA=1;            end;
  if ~isfield(options,'Bandwidth')      
    options.Bandwidth=2000;      end;
  if ~isfield(options,'FieldStrength')  
    options.FieldStrength=3;     end;
  if ~isfield(options,'Matrix')         
    options.Matrix=64;           end;
  if ~isfield(options,'SliceThickness') 
    options.SliceThickness=3;    end;
  if ~isfield(options,'TE')             
    options.TE=30;               end;

  if ~isfield(options,'AddNoise')       
    options.AddNoise=true;       end;
  if ~isfield(options,'FatSat')         
    options.FatSat=false;        end;
  if ~isfield(options,'BlipUp')         
    options.BlipUp=false;        end;
  if ~isfield(options,'Dropout')        
    options.Dropout=true;        end;
  if ~isfield(options,'AddSinus')       
    options.AddSinus=false;      end;
  if ~isfield(options,'PhaseEncodeLR')  
    options.PhaseEncodeLR=false; end;

  % print out options
  disp(' ')
  disp(['        GRAPPA = ' sprintf('%4d',options.GRAPPA)           ...
       '       AddNoise = ' bool2str(options.AddNoise)]);
  disp(['     Bandwidth = ' sprintf('%4d',options.Bandwidth)        ...
       '         FatSat = ' bool2str(options.FatSat)]);
  disp([' FieldStrength =  ' sprintf('%2.1f',options.FieldStrength)  ...
       '         BlipUp = ' bool2str(options.BlipUp)]);
  disp(['        Matrix = ' sprintf('%4d',options.Matrix)           ...
       '        Dropout = ' bool2str(options.Dropout)]);
  disp(['SliceThickness =  ' sprintf('%2.1f',options.SliceThickness) ...
       '       AddSinus = ' bool2str(options.AddSinus)]);
  disp(['            TE = ' sprintf('%4d',options.TE)               ...
       '  PhaseEncodeLR = ' bool2str(options.PhaseEncodeLR)]);
  disp(' ')
        
  
function str = bool2str(val)
  
  if val, str='true';
  else    str='false'; end;