function [T1, T2] = get_relaxation_times(field,tissue)
%  function [T1 T2] = get_relaxation_times(field,tissue)
%
%  field can be 1 1.5, 3, 4 or 7; tissue can be 'wm', 'gm', or 'blood'
%  Original author: Karla M
%  Last update : Olivia V, 24.03.2016
%  CAUTION: Literature values can vary strongly
%
%  ADDING INSTRUCTIONS:
%  Apend the relevant T1/T2 values under the field strength case and apend the
%  first author name and year in the comment. Make sure the order of values
%  and references matches. Put your name in the last update.

% now specify defaults before checking for other params
    switch (field)     
     case {1.5}
       T1gm = [900 1220 1197 1188];  % T1s from [breger 1989, Yong-Hing 2005,Wright 2008, Rooney 2007]   
       T2gm = [110 80];              % T2s [Scaled from wansapura 1999, Yong-Hing 2008]
       T1wm = [650 646 840 681];     % T1s from [breger 1989, Wright 2008, Yong-Hing 2005, Schmierer 2008];   
       T2wm = [90 80 80 81];         % T2s from [breger 1989; Cox 2010, Yong-Hing 2005, Schmierer 2008]
       T1b  = [1540];                % T1s from [Rooney 2007]
       T2b  = [250];                 % T2s from ? 
     case {3}
       T1gm = [1330 1763 1607];      % T1s from [Wansapura 1999, Gelman 2001, Wright 2008
       T2gm = [90 76];               % T2s from [Wansapura 1999, Cox 2010]
       T1wm = [830 847 838] ;        % T1s from [Wansapura 1999, Gelman 2001, Wright 2008]
       T2wm = [80 71];               % T2s from [Wansapura 1999,Cox 2010
       T1b  = [1400];                % T1s from ?
       T2b  = [160];                 % T2s from ?
     case {7}
       T1gm = [1650 2132 1939];      % T1s from [Wright 2006, Rooney 2007, Wright 2008]   
       T2gm = [60 47 46];            % T2s from [Pfeuffer 2004, Cox 2010 (frontal and occipital)] 
       T1wm = [1020 1220 1126];      % T1s from [Wright 2006,Rooney 2007, Wright 2008] 
       T2wm = [55 47];               % T2s from [Pfeuffer 2004, Cox 2010]
       T1b  = [2587];                % T1s from [Rooney 2007]
       T2b  = [30];                  % T2s from ?
      otherwise
       error('unknown field strength specified');
    end
%      ?
    switch (tissue)
     case{'gm'}
       T1=mean(T1gm); T2=mean(T2gm);
     case{'wm'}
       T1=mean(T1wm); T2=mean(T2wm);
     case{'blood'}
       T1=mean(T1b); T2=mean(T2b);
     otherwise
       error('unknown tissue type specified');
    end
end
