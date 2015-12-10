function varargout = dMRIsim(varargin)
% DTI tutorial
% wenchuan@fmrib, yuhang@fmrib

% DMRISIM MATLAB code for dMRIsim.fig
%      DMRISIM, by itself, creates a new DMRISIM or raises the existing
%      singleton*.
%
%      H = DMRISIM returns the handle to a new DMRISIM or the handle to
%      the existing singleton*.
%
%      DMRISIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DMRISIM.M with the given input arguments.
%
%      DMRISIM('Property','Value',...) creates a new DMRISIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dMRIsim_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dMRIsim_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dMRIsim

% Last Modified by GUIDE v2.5 20-Nov-2014 21:16:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dMRIsim_OpeningFcn, ...
                   'gui_OutputFcn',  @dMRIsim_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before dMRIsim is made visible.
function dMRIsim_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dMRIsim (see VARARGIN)

% Choose default command line output for dMRIsim
handles.output = hObject;

% wenchuan@fmrib : initialse parameters
handles.TE = 90;
% handles.thick = 1.25;
% handles.res = 1.25;

handles.thick = 2;
handles.res = 2;




handles.field = 1.5;
handles.BW = 1000;
handles.bvals = 1000;
handles.nodiffdir = 7;
handles.prepmethod = 'No_eddy_currents';
handles.phaseDir = 'AP';
handles.blipDir = 'Down';
handles.imshowcontrast = 1;

handles.DataReady = 0;
handles.dir_changed = 0;
handles.scrollworked = 0;
handles.voxelselworked = 0;
handles.contrastscrollworked = 0;
handles.intensity_scale_max = 1;
handles.b0indx = 1;

handles.TR = 10;

% load('b0_img');
% load('dwiimg');

% handles.b0img = b0_img;
% handles.dwiimg = dwiimg;

% setenv('PATH','/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/fsl/bin')
% 
% % set FSL environment
% setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is
% setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be

% loadin bvecs and bvals
% [fid,errstr]=fopen('bvals');
% tline = fgetl(fid);
% bvals_all = str2num(tline);
% fclose(fid);
% handles.nodiffdir = size(bvals_all,2);
% 
% bvecs_all = zeros(3,handles.nodiffdir);
% 
% [fid,errstr]=fopen('bvecs');
% for ii = 1 : 3
%     tline = fgetl(fid);
%     bvecs_all(ii,:) = str2num(tline);
% end
% fclose(fid);
% handles.bvecs_all = bvecs_all;
% handles.bvals_all = bvals_all;

addpath NIfTI

axes(handles.axes1);
cla
axis off

axes(handles.axes2);
cla
axis off

axes(handles.axes3);
cla
axis off

axes(handles.axes4);
cla
axis off

axes(handles.axes8);
cla
axis off


axes(handles.axes9);
cla
axis off

set(handles.voxel_sel,'Enable','off');
set(handles.DiffVol_tag,'Enable','off');
set(handles.Contrast_tag,'Enable','off');


handles.x = 1;
handles.y = 1;



% handles.axes1 = imshow(abs(handles.b0img),[]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dMRIsim wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dMRIsim_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;







% TE is affected by b value and diffusion preparation method

function TE_time_Callback(hObject, eventdata, handles)
% hObject    handle to TE_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TE_time as text
%        str2double(get(hObject,'String')) returns contents of TE_time as a double

TE_val = str2double(get(hObject, 'String'));
TE_val=round(TE_val);

% minTE=round(8990*(handles.BW^-0.78)); % ?? why

% very simple simulation
if (strcmp(handles.prepmethod,'No_eddy_currents'))||(strcmp(handles.prepmethod,'Monopolar'))

   minTE = handles.bvals * (10/1000) + 70;
    
elseif (strcmp(handles.prepmethod,'Bipolar'))
   
   minTE = handles.bvals * (10/1000) + 70;
   minTE = 1.1 * minTE; 
end
    
% minTE = 70; % if no sense applied

%check that it is a number
if isnan(TE_val)
    TE_val=handles.TE; %if not a valid entry keep old value
end
%ensure TE is between minimum and 200 ms
if TE_val < minTE
    TE_val=minTE;
end
if TE_val >300
    TE_val=300;
end



%convert to string for display
val=int2str(TE_val);
%update TE in UI with new value (if it has changed)
set(hObject, 'String', val);

%store TE to handle
handles.TE=TE_val;
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function TE_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TE_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function Res_Callback(hObject, eventdata, handles)
% hObject    handle to Res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Res as text
%        str2double(get(hObject,'String')) returns contents of Res as a double

%germuska@fmrib: read Resolution value from UI
Res_val = str2double(get(hObject, 'String'));

%check that it is a number
if isnan(Res_val)
    Res_val=2;
end
%ensure Resolution is between 1.25 and 20 mm
if Res_val < 1.25
    Res_val=1.25;
end
if Res_val >10
    Res_val=10;
end
%convert to string for display
val=num2str(Res_val);
%update TE in UI with new value (if it has changed)
set(hObject, 'String', val);

%store Resolution to handle
handles.res=Res_val;
% Update handles structure
guidata(hObject, handles);






% --- Executes during object creation, after setting all properties.
function Res_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SlcThk_Callback(hObject, eventdata, handles)
% hObject    handle to SlcThk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SlcThk as text
%        str2double(get(hObject,'String')) returns contents of SlcThk as a double
Thick_val = str2double(get(hObject, 'String'));

%check that it is a number
if isnan(Thick_val)
    Thick_val=2;
end
%ensure Thickness is between 0.5 and 20 mm
if Thick_val < 1.25
    Thick_val=1.25;
end
if Thick_val >10
    Thick_val=10;
end
%convert to string for display
val=num2str(Thick_val);
%update TE in UI with new value (if it has changed)
set(hObject, 'String', val);

%store Thickness to handle
handles.thick=Thick_val;
% Update handles structure
guidata(hObject, handles);





% --- Executes during object creation, after setting all properties.
function SlcThk_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SlcThk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fieldstr_Callback(hObject, eventdata, handles)
% hObject    handle to fieldstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fieldstr as text
%        str2double(get(hObject,'String')) returns contents of fieldstr as a double


val = get(hObject, 'Value');
switch val;
    case 1
        handles.field=1.5;
    case 2 
        handles.field=3;
    case 3
        handles.field=7;
end

% Update handles structure
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function fieldstr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bw_pixel_Callback(hObject, eventdata, handles)
% hObject    handle to bw_pixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bw_pixel as text
%        str2double(get(hObject,'String')) returns contents of bw_pixel as a double

%germuska@fmrib: read TE value from UI
BW_val = str2double(get(hObject, 'String'));
%round to an int
BW_val=round(BW_val);
%check that it is a number
if isnan(BW_val)
    BW_val=handles.BW; %keep old bandwidth if not a valid number
end
%ensure BW is between 750 and 2500 ms %germuska@fmrib: 2012 changed limits
%to be more realistic

if BW_val < 750
    BW_val=750;
end
if BW_val >2500
    BW_val=2500;
end


%convert to string for display
val=int2str(BW_val);
%update BW in UI with new value (if it has changed)
set(hObject, 'String', val);



%germuska@fmrib: edit 2012: update TE value to make sure it is still valid

TE_val=handles.TE;
minTE=round(8990*(BW_val^-0.78));
if TE_val<minTE
    stringTE=num2str(minTE); 
    set(handles.TE_time,'String',stringTE);
end

%store BW to handle
handles.BW=BW_val;
% Update handles structure
guidata(hObject, handles);





% --- Executes during object creation, after setting all properties.
function bw_pixel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bw_pixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bvals_Callback(hObject, eventdata, handles)
% hObject    handle to bvals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bvals as text
%        str2double(get(hObject,'String')) returns contents of bvals as a double

bvals = str2double(get(hObject, 'String'));
bvals = round(bvals);
%check that it is a number
if isnan(bvals)
    bvals=handles.bvals; %if not a valid entry keep old value
end
%ensure TE is between minimum and 200 ms
if bvals < 0
    bvals=0;
end
if bvals >20000
    bvals=20000;
end



%convert to string for display
val=int2str(bvals);
%update in UI with new value (if it has changed)
set(hObject, 'String', val);

%store to handle
handles.bvals=bvals;

% update TE

if (strcmp(handles.prepmethod,'No_eddy_currents'))||(strcmp(handles.prepmethod,'Monopolar'))
   minTE = handles.bvals * (10/1000) + 70;    
elseif (strcmp(handles.prepmethod,'Bipolar'))   
   minTE = handles.bvals * (10/1000) + 70;
   minTE = 1.1 * minTE; 
end

TE_val = get(handles.TE_time, 'String');
TE_val = str2num(TE_val);
if TE_val < minTE
    TE_val=minTE;
end
if TE_val >300
    TE_val=300;
end

%convert to string for display
val=int2str(TE_val);
%update TE in UI with new value (if it has changed)
set(handles.TE_time, 'String', val);

%store TE to handle
handles.TE=TE_val;


% Update handles structure
guidata(hObject, handles);





% --- Executes during object creation, after setting all properties.
function bvals_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bvals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in no_diffdir.
function no_diffdir_Callback(hObject, eventdata, handles)
% hObject    handle to no_diffdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns no_diffdir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from no_diffdir


val = get(hObject, 'Value');

% val = str2num(val);

% no_dir_vector = [7,12,21,22,25,32,33,35,64,68,136,139,140,142];

no_dir_vector = [7,21,22,25,32,33,35,64,68,136,139,140,142];

nodir = no_dir_vector(val);


handles.dir_changed = 1;


%store Resolution to handle
handles.nodiffdir=nodir;
% Update handles structure
guidata(hObject, handles);






% --- Executes during object creation, after setting all properties.
function no_diffdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to no_diffdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in diffprep.
function diffprep_Callback(hObject, eventdata, handles)
% hObject    handle to diffprep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns diffprep contents as cell array
%        contents{get(hObject,'Value')} returns selected item from diffprep

val = get(hObject, 'Value');
switch val;
    
    case 1
        handles.prepmethod = 'No_eddy_currents';
    case 2
        handles.prepmethod = 'Monopolar';
    case 3
    
        handles.prepmethod = 'Bipolar';
end



if (strcmp(handles.prepmethod,'No_eddy_currents'))||(strcmp(handles.prepmethod,'Monopolar'))
   minTE = handles.bvals * (10/1000) + 70;    
elseif (strcmp(handles.prepmethod,'Bipolar'))   
   minTE = handles.bvals * (10/1000) + 70;
   minTE = 1.1 * minTE; 
end

TE_val = get(handles.TE_time, 'String');
TE_val = str2num(TE_val);

if TE_val < minTE
    TE_val=minTE;
end
if TE_val >300
    TE_val=300;
end

%convert to string for display
val=int2str(TE_val);
%update TE in UI with new value (if it has changed)
set(handles.TE_time, 'String', val);

%store TE to handle
handles.TE=TE_val;









guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function diffprep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diffprep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end













% --- Executes on key press with focus on scanbutton and none of its controls.
function scanbutton_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to scanbutton (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


key = get(gcf,'CurrentKey');
    if(strcmp (key , 'return'))
%         pushbutton1_Callback(hObject, eventdata, handles)
        scanbutton_Callback(hObject, eventdata, handles)
        guidata(hObject, handles);
    end





% --- Executes on button press in scanbutton.
function scanbutton_Callback(hObject, eventdata, handles)
% hObject    handle to scanbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


disp('Scanning, please wait');
[bvalall,diffdirall] = generate_bval_dir(handles);
disp('generate diffusion weighted images');
handles.bvecs_all = diffdirall;
handles.bvals_all = bvalall;
tmpdir = 'mydatadir'; 

res = load_bedpostx_results(tmpdir);
sig = basgen(bvalall,diffdirall,res.s0,res.d,res.f,res.th,res.ph);

load fmap;

% 1.25*1.25*1.25    145*174*10
[nx1,ny1,nz1] = size(res.mask);
sig = reshape(sig,nx1,ny1,nz1,handles.nodiffdir);

% load the partial volume map

p_csf = load_nii([tmpdir '/T1w_acpc_dc_restore_1_pve_0.nii.gz']);
p_csf = p_csf.img;

p_gm = load_nii([tmpdir '/T1w_acpc_dc_restore_1_pve_1.nii.gz']);
p_gm = p_gm.img;

p_wm = load_nii([tmpdir '/T1w_acpc_dc_restore_1_pve_2.nii.gz']);
p_wm = p_wm.img;

% t2 value for wm, gm , csf. based on 3T
t2_wm = 69;
t2_gm = 100;
t2_csf = 2000;

te_ref = 89;

decay_csf = p_csf*exp(-(handles.TE - te_ref)/t2_csf);
decay_csf = repmat(decay_csf,[1,1,1,handles.nodiffdir]);

decay_wm = p_wm*exp(-(handles.TE - te_ref)/t2_wm);
decay_wm = repmat(decay_wm,[1,1,1,handles.nodiffdir]);

decay_gm = p_gm*exp(-(handles.TE - te_ref)/t2_gm);
decay_gm = repmat(decay_gm,[1,1,1,handles.nodiffdir]);


sig = sig.*(decay_csf + decay_wm + decay_gm);

disp('image resampling');

xo = (1:nx1)*1.25;
yo = (1:ny1)*1.25;
zo = (1:nz1)*1.25;

[xo,yo,zo] = meshgrid(yo,xo,zo);

xn = handles.res : handles.res : nx1*1.25;
yn = handles.res : handles.res : ny1*1.25;
zn = handles.thick : handles.thick : nz1*1.25;

[xn,yn,zn] = meshgrid(yn,xn,zn);
signew = zeros(size(xn));
signew = repmat(signew,[1,1,1,handles.nodiffdir]);
for i = 1 : handles.nodiffdir
signew(:,:,:,i)   = interp3(xo,yo,zo,sig(:,:,:,i),xn,yn,zn);
end


slabcenter = ceil(size(xn,3)/2);
data = squeeze(signew(:,:,slabcenter,:));

mskdata = zeros(size(data(:,:,1)));
mskdata(find(data(:,:,1))) = 1;

% interpolate fmap

xo = (1:nx1)*1.25;
yo = (1:ny1)*1.25;
[xo,yo] = meshgrid(yo,xo);
xn = handles.res : handles.res : nx1*1.25;
yn = handles.res : handles.res : ny1*1.25;

[xn,yn] = meshgrid(yn,xn);
fmapn = interp2(xo,yo,fmap,xn,yn);
fmapn = fmapn.*mskdata;
% fmapn = fmapn/2/pi;
% for B0 induced distortion, it is consistant across all b value all
% directions
MaxUpSampFactor = 15;
UpSampFactor = round(min(MaxUpSampFactor,21/handles.res));

% scaled according to field strength;
fmapn = fmapn/3 * handles.field;


% merge the field map with eddy field map/ but when this is no eddy
% simulation, still run this

if (strcmp(handles.prepmethod,'No_eddy_currents') == 1)
 data = Distortion_dwiEPI(data, fmapn, UpSampFactor, handles.TE, handles.BW, size(fmapn), handles.phaseDir, handles.blipDir);
end

disp('signal scaling and adding noise');

% change the signal intensity according to res, slicethickness, bw,
% fieldstr, te(partial volume model)

res_ref = 1.25;
slithick_ref = 1.25;
fieLdstr_ref = 3;
bw_ref = 1490;
readout_length = 174;


% find a b0 image
b0ind = find(handles.bvals_all == 0);
b0imgtmp = data(:,:,b0ind(1));
handles.b0indx = b0ind(1);
%%
% Add thermal noise, scaled by bandwidth
ThermalNoiseSD = (size(b0imgtmp,2) * handles.BW) / (readout_length * bw_ref) * 600;

% Add physiological noise - cardiac plus random
%?? spatial dependence ??

MeanSig = mean(abs(b0imgtmp(:)));
FieldStrengths = [1.5 3 7];
lambdas = [0.0123 0.0107 0.0086]; % From Triantafyllou, Neuroimage 2005
lambda = interp1(FieldStrengths,lambdas,handles.field);
PhysiolNoiseSD = lambda * MeanSig * 60;

%%

added_noise = PhysiolNoiseSD*rand(size(b0imgtmp)) + ThermalNoiseSD * rand(size(b0imgtmp)); 



signal_scaling_factor = (handles.res/res_ref)^2 * handles.thick/slithick_ref...
    *handles.field/fieLdstr_ref;

data = data * signal_scaling_factor;
data_freenoise = data;
% data = data + repmat(added_noise,[1,1,handles.nodiffdir]);


if handles.DataReady==0
    pause(3);
end



handles.dwiimg = data;

% distorted the images by eddy current

%  eddy_phase_map_all = zeros(size(data));

if (strcmp(handles.prepmethod,'No_eddy_currents') == 0)

    
    disp('eddy current simulation start');
    
   referx = @(xx,timedelay)((0.00188006*exp(-timedelay/0.883301) ...
    +0.00394905*exp(-timedelay/0.202146) ...
    +0.000998472*exp(-timedelay/0.00885462) ...
    -0.00159343*exp(-timedelay/0.00572957)...
    +0.00151317*exp(-timedelay/0.00199695))*xx);
 
   refery= @(yy,timedelay) ((0.00203354*exp(-timedelay/1.40515) ...
            +0.00782493*exp(-timedelay/0.383857) ...
            +0.00259124*exp(-timedelay/0.102849) ...
            +0.00543571*exp(-timedelay/0.00168112)...
            -0.00446654*exp(-timedelay/0.00121779))*yy);
      

   ksp_new = zeros(size(data));
   
  
        
end


if (strcmp(handles.prepmethod,'Monopolar') == 1)
  
    TE_val = handles.TE;
    
    time_decay1 = (TE_val - (TE_val/2 - 18))/1000;
    time_decay2 = (TE_val - (TE_val/2 -7))/1000;
    time_decay3 = (TE_val - (TE_val/2+8))/1000;
    time_decay4 = (TE_val - (TE_val/2 + 19))/1000;
    
 
    
    
% tic
for idir = 1 : handles.nodiffdir
    
    
    disp(['eddy current simulation for direction ', num2str(idir)]);
    
    b_tmp = handles.bvals_all(idir);
    vec_tmp = handles.bvecs_all(:,idir);
    
    g_amp_tmp = 10/sqrt(1000/b_tmp);
    
    g_slr_tmp_x = g_amp_tmp * vec_tmp(1)/100;
    g_slr_tmp_y = g_amp_tmp * vec_tmp(2)/100;
    
    eddy_x = -1 * referx(g_slr_tmp_x,time_decay1) + referx(g_slr_tmp_x,time_decay2) + ...
        -1 * referx(g_slr_tmp_x,time_decay3) + referx(g_slr_tmp_x,time_decay4);
    
    eddy_y = -1 * refery(g_slr_tmp_y,time_decay1) + refery(g_slr_tmp_y,time_decay2) + ...
        -1 * refery(g_slr_tmp_y,time_decay3) + refery(g_slr_tmp_y,time_decay4);
    
    
    % eddy is too strong..
    eddy_x = eddy_x/10;
    eddy_y = eddy_y/10;
    
%     eddy_x = 0;
%     eddy_y = 0;
    



     %% wenchuan new trying 
    [xsz,ysz] = size(squeeze(data(:,:,1)));
     eddy_phase_map = zeros(xsz,ysz);
     
     xcenter = floor(xsz/2);
     ycenter = floor(ysz/2);
     xind = ((1:xsz) - xcenter)*handles.res/1000;
     yind = ((1:ysz) - ycenter)*handles.res/1000;
     
     for iix = 1 : xsz
         for iiy = 1 : ysz
             
           eddy_phase_map(iix,iiy) = 42.576*10^6 * (eddy_x * xind(iix)  + eddy_y * yind(iiy));
         end
     end
     
     data(:,:,idir) = Distortion_dwiEPI(data(:,:,idir), (eddy_phase_map + fmapn), UpSampFactor, handles.TE, handles.BW, size(eddy_phase_map), handles.phaseDir, handles.blipDir);
   
    
    
    
    
end
  
 
% toc
% handles.dwiimg = ifft2c(ksp_new); 

handles.dwiimg = data; 
    
    
elseif strcmp(handles.prepmethod, 'Bipolar') == 1


 
  TE_val = handles.TE;
    
    time_decay1 = (TE_val - (TE_val/4 - 13))/1000;
    time_decay2 = (TE_val - (TE_val/4 -7))/1000;
    time_decay3 = (TE_val - (TE_val/4 +8))/1000;
    time_decay4 = (TE_val - (TE_val/4 + 19))/1000;  
    time_decay5 = (TE_val - (3*TE_val/4 - 18))/1000;  
    time_decay6 = (TE_val - (3*TE_val/4 - 7))/1000;  
    time_decay7 = (TE_val - (3*TE_val/4 + 8))/1000;  
    time_decay8 = (TE_val - (3*TE_val/4 + 14))/1000;  
    
    
      
for idir = 1 : handles.nodiffdir
    
    disp(['eddy current simulation for direction ', num2str(idir)]);
    
    b_tmp = handles.bvals_all(idir);
    vec_tmp = handles.bvecs_all(:,idir);
    
    g_amp_tmp = 10/sqrt(1000/b_tmp);
    
   
    g_slr_tmp_x = g_amp_tmp * vec_tmp(1)/100;
    g_slr_tmp_y = g_amp_tmp * vec_tmp(2)/100;
    
    % sign of eddy : - + + - - + + -
    
    eddy_x = -1 * referx(g_slr_tmp_x,time_decay1) + referx(g_slr_tmp_x,time_decay2) + ...
          referx(g_slr_tmp_x,time_decay3) + (-1) * referx(g_slr_tmp_x,time_decay4) + ...
        + (-1) * referx(g_slr_tmp_x,time_decay5) + referx(g_slr_tmp_x,time_decay6) + ...
        referx(g_slr_tmp_x,time_decay7) + (-1) * referx(g_slr_tmp_x,time_decay8);
    
    eddy_y = -1 * refery(g_slr_tmp_y,time_decay1) + refery(g_slr_tmp_y,time_decay2) + ...
         refery(g_slr_tmp_y,time_decay3) + (-1) * refery(g_slr_tmp_y,time_decay4) + ... 
        (-1) * refery(g_slr_tmp_y,time_decay5) + refery(g_slr_tmp_y,time_decay6) + ...
         refery(g_slr_tmp_y,time_decay7) + (-1) * refery(g_slr_tmp_y,time_decay8);
    

    eddy_x = eddy_x/10;
    eddy_y = eddy_y/10;
    




     %% wenchuan new trying 
    [xsz,ysz] = size(squeeze(data(:,:,1)));
     eddy_phase_map = zeros(xsz,ysz);
     
     xcenter = floor(xsz/2);
     ycenter = floor(ysz/2);
     xind = ((1:xsz) - xcenter)*handles.res/1000;
     yind = ((1:ysz) - ycenter)*handles.res/1000;
     
     for iix = 1 : xsz
         for iiy = 1 : ysz
             
           eddy_phase_map(iix,iiy) = 42.576*10^6 * (eddy_x * xind(iix)  + eddy_y * yind(iiy));
         end
     end
     
     data(:,:,idir) = Distortion_dwiEPI(data(:,:,idir), (eddy_phase_map + fmapn), UpSampFactor, handles.TE, handles.BW, size(eddy_phase_map), handles.phaseDir, handles.blipDir);

    
 
    
 end
      
    
    handles.dwiimg = data; 
    
end
handles.dwiimg = handles.dwiimg + repmat(added_noise,[1,1,handles.nodiffdir]);
% scale the dwiimg to show it properly 
handles.intensity_scale_max = max(abs(handles.dwiimg(:)))/2;





disp('dti fitting');



% res2 = fast_dtifit_eig(data,bvalall,diffdirall,data_freenoise);

% data_for_mask = handles.dwiimg(:,:,b0ind(1));
data_freenoise = ones(size(data_freenoise));
res2 = fast_dtifit_eig(abs(handles.dwiimg),bvalall,diffdirall,data_freenoise);


Mean_adc = res2.MD;

dti_V1 = res2.V1;
dti_FA = res2.FA;

V1_mod_FA = abs(dti_V1 .* repmat(dti_FA,[1,1,3]));
V1_mod_FA = V1_mod_FA/max(V1_mod_FA(:));
V1_mod_FA_forshow = zeros(size(V1_mod_FA,2),size(V1_mod_FA,1),3);
for ii = 1 : 3
    V1_mod_FA_forshow(:,:,ii) = rot90(V1_mod_FA(:,:,ii));
end

handles.DataReady = 1;

 no_of_volumes = handles.nodiffdir;
 set(handles.DiffVol_tag,'Max',no_of_volumes-1);
set(handles.DiffVol_tag,'Min',0);

set(handles.DiffVol_tag, 'SliderStep', [1/(no_of_volumes-1) , 4/(no_of_volumes-1) ]);


if handles.dir_changed == 1 
   handles.scrollworked = 0;
   set(handles.DiffVol_tag, 'Value', 0); 
   handles.dir_changed = 0;
end



guidata(hObject, handles);





axes(handles.axes1);



cla
%imshow(abs(rot90(handles.dwiimg(:,:,b0ind(1)))),[0 handles.intensity_scale_max * handles.imshowcontrast]);axis image; axis off;
imagesc(abs(rot90(handles.dwiimg(:,:,b0ind(1)))),[0 handles.intensity_scale_max * handles.imshowcontrast]);axis image; axis off;

axes(handles.axes2);
cla
if handles.scrollworked == 0
%imshow(abs(rot90(handles.dwiimg(:,:,1))),[0 handles.intensity_scale_max * handles.imshowcontrast]);axis image; axis off;
imagesc(abs(rot90(handles.dwiimg(:,:,1))),[0 handles.intensity_scale_max * handles.imshowcontrast]);axis image; axis off; % RF added
else    
%imshow(abs(rot90(handles.dwiimg(:,:,handles.diffvol))),[0 handles.intensity_scale_max * handles.imshowcontrast]);axis image; axis off;        
imagesc(abs(rot90(handles.dwiimg(:,:,handles.diffvol))),[0 handles.intensity_scale_max * handles.imshowcontrast]);axis image; axis off; % RF added     
end

if handles.voxelselworked == 1
    
rectangle('Position',[(handles.x)-2 (handles.y)-2  5 5], 'LineWidth',2, 'EdgeColor','r');
    
end


if (handles.voxelselworked == 1)
[~,n2,~] = size(handles.dwiimg);    
% the x, y is reordered
sigforshow = squeeze(abs(handles.dwiimg(handles.x,(n2-handles.y+1),:)));
sigforshow = sigforshow/sigforshow(1);

axes(handles.axes9);
cla
plot(sigforshow,'.');

if handles.scrollworked == 1
hold on
plot(handles.diffvol,sigforshow(handles.diffvol),'r*'); 
end
end






axes(handles.axes3);
cla
%imshow(abs(rot90(Mean_adc)),[0 0.002]);axis image; axis off;
imagesc(abs(rot90(Mean_adc)),[0 0.002]);axis image; axis off; % RF added
colormap('gray') % RF added

axes(handles.axes4);
cla
% imshow(rot90(fa),[]);axis image; axis off;
image((V1_mod_FA_forshow));axis image; axis off;


% show the diffusion direction for the first volume
axes(handles.axes8);
cla
stvec = [2 1 0;1 2 0; 1 1 -1];
endvec = [0 1 0; 1 0 0; 1 1 1];

set(gca,'ColorOrder',abs(endvec - stvec)/2); % divide 2 to make it in [0 1]

% arrow3(stvec,endvec,'k',1,3,0.1);

arrow3(stvec,endvec,'o',1,3,0.1);

text(-0.2,1,-0.1,'Right');
text(2.2,1,-0.1,'Left');
text(0.8,2.6,-0.1,'Posterior');
text(0.8,-0.2,0,'Anterior');
text(0.8,1,-1.2,'Inferior');
text(0.8,1,1.2,'Superior');
 
 hold on
 
 if handles.scrollworked == 0
vec2s = (diffdirall(:,1)).';
 else
   vec2s = (diffdirall(:,handles.diffvol)).';  
 end
 set(gca,'ColorOrder',abs(vec2s))
 
 vec2s(1) = -vec2s(1);
 vec2s(2) = -vec2s(2);
 vec2s = vec2s + [1 1 0];
 if(vec2s ~= [1 1 0]);


arrow3([1 1 0],vec2s,'o',3,4);

 end
axis off

disp('Done!');


% update the dwi image info

if handles.scrollworked == 0

outstring = ['b = ' num2str(handles.bvals_all(1)) ' s/mm2'];

set(handles.bvals_show,'String',outstring);

outstring = ['direction = ' num2str(handles.bvecs_all(1,1)) '  ' num2str(handles.bvecs_all(2,1)) '  ' num2str(handles.bvecs_all(3,1))];


set(handles.bvecs_show,'String',outstring);

else
    
outstring = ['b = ' num2str(handles.bvals_all(handles.diffvol)) ' s/mm2'];

set(handles.bvals_show,'String',outstring);

outstring = ['direction = ' num2str(handles.bvecs_all(1,handles.diffvol)) '  ' num2str(handles.bvecs_all(2,handles.diffvol)) '  ' num2str(handles.bvecs_all(3,handles.diffvol))];

set(handles.bvecs_show,'String',outstring);
    
end


% show scan time
scantime = handles.TR * handles.nodiffdir;
scantime_min = floor(scantime/60);
scantime_sec = scantime - scantime_min*60;
 stringscantim=[num2str(scantime_min),' : ',num2str(scantime_sec)]; 
  set(handles.scantime,'String',stringscantim);



set(handles.voxel_sel,'Enable','on');
set(handles.DiffVol_tag,'Enable','on');
set(handles.Contrast_tag,'Enable','on');








% --- Executes on slider movement.
function DiffVol_tag_Callback(hObject, eventdata, handles)
% hObject    handle to DiffVol_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.scrollworked = 1; % before the first try


handles.diffvol = get(hObject,'Value');

% handles.diffvol = round(10*handles.diffvol);
handles.diffvol = round(handles.diffvol + 1);




outstring = num2str(handles.diffvol);

set(handles.Vol_index,'String',outstring);

if handles.DataReady > 0
    
axes(handles.axes2);
cla
%imshow(abs(rot90(handles.dwiimg(:,:,handles.diffvol))),[0 handles.intensity_scale_max * handles.imshowcontrast]);axis image; axis off; 
imagesc(abs(rot90(handles.dwiimg(:,:,handles.diffvol))),[0 handles.intensity_scale_max * handles.imshowcontrast]);axis image; axis off; % RF added
colormap('gray') % RF added

if handles.voxelselworked == 1

% imshow(im_one_voxel,[]);
rectangle('Position',[(handles.x)-2 (handles.y)-2  5 5], 'LineWidth',2, 'EdgeColor','r');

end
% update the plot of signal



if (handles.voxelselworked == 1)
[~,n2,~] = size(handles.dwiimg);    
% the x, y is reordered
sigforshow = squeeze(abs(handles.dwiimg(handles.x,(n2-handles.y+1),:)));
sigforshow = sigforshow/sigforshow(1);


axes(handles.axes9);
cla
plot(sigforshow,'.');

hold on
plot(handles.diffvol,sigforshow(handles.diffvol),'r*');         
end



axes(handles.axes8);
cla
stvec = [2 1 0;1 2 0; 1 1 -1];
endvec = [0 1 0; 1 0 0; 1 1 1];

set(gca,'ColorOrder',abs(endvec - stvec)/2);

%  arrow3(stvec,endvec,'k',1,3,0.1);

 arrow3(stvec,endvec,'o',1,3,0.1);
 
text(-0.2,1,-0.1,'Right');
text(2.2,1,-0.1,'Left');
text(0.8,2.6,-0.1,'Posterior');
text(0.8,-0.2,0,'Anterior');
text(0.8,1,-1.2,'Inferior');
text(0.8,1,1.2,'Superior');
 
 
 
 
 hold on
vec2s = (handles.bvecs_all(:,handles.diffvol)).';
set(gca,'ColorOrder',abs(vec2s))


 vec2s(1) = -vec2s(1);
 vec2s(2) = -vec2s(2);
 vec2s = vec2s + [1 1 0];
 
 vec1 = [1 1 0];
 if(sum(abs(vec2s - vec1)) ~= 0)
% arrow3(vec1,vec2s,'b',3,4);

% arrow3(zeros(4,2),[10*rand(4,1),500*rand(4,1)],'o*/',w,h,0)
    
arrow3([1 1 0],vec2s,'o',3,4);
 end
axis off







% 
outstring = ['b =  ' num2str(handles.bvals_all(handles.diffvol)) '  s/mm2'];
% 
set(handles.bvals_show,'String',outstring);

outstring = ['direction =  [' num2str(handles.bvecs_all(1,handles.diffvol))...
    '   ' num2str(handles.bvecs_all(2,handles.diffvol)) '   ' num2str(handles.bvecs_all(3,handles.diffvol)) ']'];

set(handles.bvecs_show,'String',outstring);


    
end



guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function DiffVol_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DiffVol_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in voxel_sel.
function voxel_sel_Callback(hObject, eventdata, handles)
% hObject    handle to voxel_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.voxelselworked = 1;

axes(handles.axes2);   % plot the place
% rectangle('visible','off');
cla


if handles.scrollworked == 1 
%imshow(abs(rot90(handles.dwiimg(:,:,handles.diffvol))),[0 handles.intensity_scale_max .* handles.imshowcontrast]);axis image; axis off; 
imagesc(abs(rot90(handles.dwiimg(:,:,handles.diffvol))),[0 handles.intensity_scale_max .* handles.imshowcontrast]);axis image; axis off; % RF added 
else
%imshow(abs(rot90(handles.dwiimg(:,:,handles.b0indx))),[0 handles.intensity_scale_max .* handles.imshowcontrast]);axis image; axis off; 
imagesc(abs(rot90(handles.dwiimg(:,:,handles.b0indx))),[0 handles.intensity_scale_max .* handles.imshowcontrast]);axis image; axis off; % RF added
end    



axes(handles.axes2);
[y,x] = ginput(1);

x=round(x);
y=round(y);

[n1,n2,nd] = size(handles.dwiimg);


if x < 1
    x = 1;
end
if y < 1
    y = 1;
end
if x > n2
    x = n2;
end
if y > n1
    y = n1;
end

sigforshow = squeeze(abs(handles.dwiimg(y,(n2-x+1),:)));
sigforshow = sigforshow/sigforshow(1);

axes(handles.axes9);
cla
plot(sigforshow,'.');


axes(handles.axes2);   % plot the place


hold on


rectangle('Position',[y-2 x-2  5 5], 'LineWidth',2, 'EdgeColor','r');

hold off
% hold on




handles.x = y;
handles.y = x;


  stringX=num2str(handles.x); 
  set(handles.x_text,'String',stringX);
  stringY=num2str(handles.y); 
  set(handles.y_text,'String',stringY);

% handles.x = x;
% handles.y = y;


guidata(hObject, handles);


% --- Executes on selection change in Phase_Dir.
function Phase_Dir_Callback(hObject, eventdata, handles)
% hObject    handle to Phase_Dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Phase_Dir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Phase_Dir

val = get(hObject, 'Value');
switch val;
    
%     handles.phaseDir = 'AP';
% handles.blipDir = 'Down';
    
    case 1
        handles.phaseDir = 'AP';
    case 2
        handles.phaseDir = 'LR';

end

guidata(hObject, handles);







% --- Executes during object creation, after setting all properties.
function Phase_Dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Phase_Dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Blip_Dir.
function Blip_Dir_Callback(hObject, eventdata, handles)
% hObject    handle to Blip_Dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Blip_Dir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Blip_Dir

val = get(hObject, 'Value');
switch val;
    
    
    case 1
        handles.blipDir = 'Down';
    case 2
        handles.blipDir = 'Up';

end

guidata(hObject, handles);





% --- Executes during object creation, after setting all properties.
function Blip_Dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Blip_Dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function Contrast_tag_Callback(hObject, eventdata, handles)
% hObject    handle to Contrast_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.contrastscrollworked = 1;

set(handles.Contrast_tag,'Max',1);
set(handles.Contrast_tag,'Min',0);

set(handles.Contrast_tag, 'SliderStep', [1/20 , 4/20 ]);

handles.imshowcontrast = get(hObject,'Value');

guidata(hObject, handles);


% change the image every time updating the scroll


axes(handles.axes1);
cla
%imshow(abs(rot90(handles.dwiimg(:,:,handles.b0indx))),[0 handles.intensity_scale_max * handles.imshowcontrast]);axis image; axis off;
imshow(abs(rot90(handles.dwiimg(:,:,handles.b0indx))),[0 handles.intensity_scale_max * handles.imshowcontrast]);axis image; axis off; % RF added

axes(handles.axes2);
cla
if handles.scrollworked == 0
imshow(abs(rot90(handles.dwiimg(:,:,1))),[0 handles.intensity_scale_max * handles.imshowcontrast]);axis image; axis off;
else    
imshow(abs(rot90(handles.dwiimg(:,:,handles.diffvol))),[0 handles.intensity_scale_max * handles.imshowcontrast]);axis image; axis off;        
end

if handles.voxelselworked == 1

% imshow(im_one_voxel,[]);
rectangle('Position',[(handles.x)-2 (handles.y)-2  5 5], 'LineWidth',2, 'EdgeColor','r');

end


% --- Executes during object creation, after setting all properties.
function Contrast_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Contrast_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
