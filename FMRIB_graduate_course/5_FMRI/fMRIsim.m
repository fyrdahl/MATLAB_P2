function varargout = fMRIsim(varargin)


%fMRI simulator written by Mike Germuska & Tom Okell. fmrib 2011.



% FMRISIM MATLAB code for fMRIsim.fig
%      FMRISIM, by itself, creates a new FMRISIM or raises the existing
%      singleton*.
%
%      H = FMRISIM returns the handle to a new FMRISIM or the handle to
%      the existing singleton*.
%
%      FMRISIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FMRISIM.M with the given input arguments.
%
%      FMRISIM('Property','Value',...) creates a new FMRISIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fMRIsim_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fMRIsim_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fMRIsim

% Last Modified by GUIDE v2.5 22-Nov-2012 09:48:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fMRIsim_OpeningFcn, ...
                   'gui_OutputFcn',  @fMRIsim_OutputFcn, ...
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


% --- Executes just before fMRIsim is made visible.
function fMRIsim_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fMRIsim (see VARARGIN)

% Choose default command line output for fMRIsim
handles.output = hObject;


%germuska@fmrib: initialise Scan parameters to default (but inappropriate) values
handles.thick=1;
handles.TE=100; %edited 2012
handles.res=1;
handles.oldres=1;
handles.field=1.5;
handles.BW=1000; %edited 2012
handles.phase='AP';
handles.blip=true;
handles.regress=false;

handles.x=1;
handles.y=1;

handles.threshold=0;
handles.DataReady=0;
handles.DataSet=0;

%Load saved data from mat files

%S=load('brain_data.mat');
S=load('sim_data_5.mat');


%store data to handles structure
handles.visual_image=S.T1_visual;
handles.visual_stats=S.zstats_1;
handles.front_image=S.T1_whole;
handles.front_stats=S.zstats_2;
handles.DropOutMap = S.DropOutMap;

handles.visual_t2star=S.T2_star_visual; %germuska@fmrib:2012 additional datasets (struct and t2 star weighted)
handles.front_t2star=S.T2_star_whole;

handles.fieldmap=S.field;


%create a graphics object containing the initial brain image
handles.brain=S.T1_visual;
%initial/dummy zstats for testing
handles.dzstats=S.zstats_1;

% Update handles structure
guidata(hObject, handles);


%germuska@fmrib: display the brain image in axes1
%handles.axes1=imshow(handles.brain);


% UIWAIT makes fMRIsim wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fMRIsim_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function TE_time_Callback(hObject, eventdata, handles)
% hObject    handle to TE_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TE_time as text
%        str2double(get(hObject,'String')) returns contents of TE_time as a double

%germuska@fmrib: read TE value from UI
TE_val = str2double(get(hObject, 'String'));
%round to an int
TE_val=round(TE_val);

%germuska@fmrib: edit for 2012: only allow valid TE values dependent on the
%bandwidth :Minimum TE approximately follows the rule 8990*Bandwith^-0.78

minTE=round(8990*(handles.BW^-0.78));

%check that it is a number
if isnan(TE_val)
    TE_val=handles.TE; %if not a valid entry keep old value
end
%ensure TE is between minimum and 200 ms
if TE_val < minTE
    TE_val=minTE;
end
if TE_val >200
    TE_val=200;
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
    Res_val=1;
end
%ensure Resolution is between 0.1 and 20 mm
if Res_val < 0.1
    Res_val=0.1;
end
if Res_val >20
    Res_val=20;
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


% --- Executes on selection change in FieldStrMenu.
function FieldStrMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FieldStrMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FieldStrMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FieldStrMenu


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
function FieldStrMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FieldStrMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SliceThick_Callback(hObject, eventdata, handles)
% hObject    handle to SliceThick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SliceThick as text
%        str2double(get(hObject,'String')) returns contents of SliceThick as a double

%germuska@fmrib: read Resolution value from UI
Thick_val = str2double(get(hObject, 'String'));

%check that it is a number
if isnan(Thick_val)
    Thick_val=1;
end
%ensure Thickness is between 0.5 and 20 mm
if Thick_val < 0.5
    Thick_val=0.5;
end
if Thick_val >20
    Thick_val=20;
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
function SliceThick_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SliceThick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ScanButton.
function ScanButton_Callback(hObject, eventdata, handles)
% hObject    handle to ScanButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



if handles.oldres~=handles.res
    handles.x=1;
    handles.y=1;
    %update display with voxel location
    stringX=num2str(handles.x); 
    set(handles.x_text,'String',stringX);
    stringY=num2str(handles.y); 
    set(handles.y_text,'String',stringY);
    handles.oldres=handles.res;
end
    
%pause to stop code from crashing on first run...
if handles.DataReady==0
    pause(3);
end


%grabData function will return z-stats
[zstats TS fit SNR tSNR MeanSig]=grabTomsData(handles);

%add zstats and timecourse data to handles structure and update gui data
handles.zstats=zstats;
handles.TS=TS;
handles.fit=fit;
handles.SNR=SNR;
handles.tSNR=tSNR;

%MeanSig=round(MeanSig./max(MeanSig(:)))*255;


func=MeanSig./max(MeanSig(:)); %germuska@fmrib 2012 add returned functional data to handles 
color_func(:,:,1)=uint8(func*255);
color_func(:,:,2)=uint8(func*255);
color_func(:,:,3)=uint8(func*255);
handles.func=color_func;


handles.DataReady=1;
% Update handles structure
guidata(hObject, handles);


%Function to upsample /(nearest neighbour) and display zstats
%as an overlay.
plotOverlay(handles)


%GRAB VOXEL DATA FROM TOM AND PLOT
grabVoxelData(handles); 

set(handles.voxel_button,'Enable','on');
set(handles.data_select,'Enable','on');
set(handles.im_weight_pop,'Enable','on');





function [zstats TS Fit SNR tSNR MeanSig]=grabTomsData(handles)


TE=handles.TE;
Resolution=handles.res;
SliceThick=handles.thick;
FieldStr=handles.field;

% Call to Tom's function 

opts.dR2star=0.0003.*handles.dzstats;

[m,n]=size(handles.dzstats);
Fmap=ones(m,n);


if handles.DataReady<1  %if first scan put appropriate t2star data in handles
    handles.t2star=handles.visual_t2star;
else
    
    val = get(handles.data_select, 'Value');
switch val;
    case 1
        handles.t2star=handles.visual_t2star;
    case 2 
        handles.t2star=handles.front_t2star;
end
    
end


%image data stored as RGB (3 channels of uint8 0 to 255)
%convert image data to numbers

%SigDensity should always be calculated from the t2star weighted data no matter what is being displayed
SigDensity=double(squeeze(handles.t2star(:,:,1))); 

%scale data to give reasonable signal density :germuska@fmrib:2012 these
%scaling values need to be updated for new input images...
opts.SigDensity=SigDensity./7;
size(opts.SigDensity);

%ones for fieldmap if visual data.%determine which field map to use + fudge
%factors to get voxel resolution looking good
if handles.DataSet==0 % Visual
    opts.FMap=Fmap;
    Resolution=Resolution/1.2;
   % opts.SigDensity=opts.SigDensity./2;
    opts.DropOutMap = [];
else % Whole brain
    opts.FMap=(handles.fieldmap)./6.28;
    opts.SigDensity=opts.SigDensity.*2;
    Resolution=Resolution/0.75;
    opts.DropOutMap = handles.DropOutMap;
end

%set matrix size from smallest dimension
if m<=n
    Smatrix=m;
else
    Smatrix=n;
end
MtxSize=round(Smatrix/Resolution);

opts.MtxSize=MtxSize;

% T.O. Set the upsampling factor used for distortion/dropout depending on the resolution
MaxUpSampFactor = 15;
opts.UpSampFactor = round(min(MaxUpSampFactor,21/Resolution));

opts.SlcThk=SliceThick;
opts.TE=TE;
opts.B0=FieldStr;
opts.BandWidthPerPix=handles.BW;
opts.PhaseEncDir= handles.phase;
opts.blip_down=handles.blip;
opts.RegressPhysioNoise=handles.regress;

[zstats SNR tSNR MeanSig TS Fit] = CalcFMRIStats(opts);

% Scale TS and Fit by MeanSig
Fit=squeeze(Fit);
TS = TS./repmat(MeanSig,[1 1 size(TS,3)])*100;
Fit = Fit./repmat(MeanSig,[1 1 size(TS,3)])*100;




%function to plot zstats overlay with specified threshold
function plotOverlay(handles)

zstats=handles.zstats;

%scale z-stats and convert to ints so that index can be converted to RGB.
zscale=255/(max(zstats(:)));

ovimg=zscale.*zstats;
ovimg=uint8(ovimg);
ovimg=double(ovimg);

%up re-sample z-stats to res of brain image

%germuska@fmrib:2012 change data set struct/functinal dependent on pop menu
%setting
val = get(handles.im_weight_pop, 'Value');
switch val;
    case 1
        [m,n,o]=size(handles.func);
        bgimg=handles.func;
    case 2 
        [m,n,o]=size(handles.brain);
        bgimg=handles.brain;
end
    

% T.O. Modified to use new function that doesn't depend on the image
% processing toolbox, Oct 2013
ovimg=toimresize(ovimg,[m n],'nearest');
zstats=toimresize(zstats,[m n],'nearest');

%convert ovimg into RGB image
ovimg=ind2rgb(ovimg,jet(255));

%create alpha transarency from zstats so that only zstats greater than
%threshold are displayed
threshold=handles.threshold;
ovimgAlphaData = zeros(m,n);
ovimgAlphaData(gt(zstats,threshold))=1;

%get axes handle and clear contents
axes(handles.axes1);
cla

%create divisions for colorbar so that it maps to zstats
%needed as zstats are now dislpayed as RGB image rather than z-stat values
a=round(100*0.9375.*max(zstats(:))/6)/100;
b=round(100*0.9375.*max(zstats(:))*2/6)/100;
c=round(100*0.9375.*max(zstats(:))*3/6)/100;
d=round(100*0.9375.*max(zstats(:))*4/6)/100;
e=round(100*0.9375.*max(zstats(:))*5/6)/100;
f=round(100*0.9375.*max(zstats(:)))/100;


%display overlay and background images
toimshow(ovimg);colorbar('YTickLabel',{a,b,c,d,e,f}); %add this line here to get colorbar to consistently display full range
hold on
toimshow(bgimg);
%hold on
I=toimshow(ovimg);
hold off
%display overimg as transparent where zstats are less than threshold
set(I,'AlphaData',ovimgAlphaData);



% --- Executes on slider movement.
function ThreshSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ThreshSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.threshold=get(hObject,'Value');

%restrict accuracy of threshold to 0.1
handles.threshold=round(10*handles.threshold)/10;
set(hObject,'Value',handles.threshold);

%update display with threhold value
outstring=num2str(handles.threshold); 
set(handles.zThresh,'String',outstring);

%if Zstats calculated update display using new threshold
if handles.DataReady>0
    plotOverlay(handles)
end

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function ThreshSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThreshSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in voxel_button.
function voxel_button_Callback(hObject, eventdata, handles)
% hObject    handle to voxel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%get input from user about voxel position


%get voxel location from user
axes(handles.axes1);
[y,x]=ginput(1);

x=round(x);
y=round(y);

%convert voxel location to zstat space. 
%germuska@fmrib 2012 determine locatino depending on displayed image
val = get(handles.im_weight_pop, 'Value');
switch val;
    case 1
        [m,n,o]=size(handles.func);
        
    case 2 
        [m,n,o]=size(handles.brain);
        
end



[p,q]=size(handles.zstats);

x=round((x/m)*p);
y=round((y/n)*q);

%check voxel location is within image bounds
if x>p
    x=p;
end
if x<1
    x=1;
end
if y>q
    y=q;
end
if y<1
    y=1;
end

handles.x=x;
handles.y=y;


%update display with voxel location
stringX=num2str(handles.x); 
set(handles.x_text,'String',stringX);
stringY=num2str(handles.y); 
set(handles.y_text,'String',stringY);

% Update handles structure
guidata(hObject, handles);


%GRAB VOXEL DATA FROM TOM AND PLOT
grabVoxelData(handles);





function grabVoxelData(handles)

x=handles.x;
y=handles.y;

%create residual data
resid_data=handles.TS-handles.fit;

%plot raw data, fits (and residual on next axes)
axes(handles.axes2);
cla
plot(squeeze(handles.TS(x,y,:)));
hold on
fit=plot(squeeze(handles.fit(x,y,:)));
set(fit,'Color','red');
xlim([0 size(handles.TS,3)])
hold off

axes(handles.axes3);
plot(squeeze(resid_data(x,y,:)));
xlim([0 size(handles.TS,3)])

%display zstat value SNR and tSNR
stringSNR=num2str(handles.SNR(x,y));
set(handles.SNR_text,'String',stringSNR);
stringtSNR=num2str(handles.tSNR(x,y));
set(handles.tSNR_text,'String',stringtSNR);
stringzstat=num2str(handles.zstats(x,y));
set(handles.zstat_text,'String',stringzstat);
    


% --- Executes on selection change in data_select.
function data_select_Callback(hObject, eventdata, handles)
% hObject    handle to data_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns data_select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from data_select


val = get(hObject, 'Value');
switch val;
    case 1
        handles.DataSet=0;
        handles.dzstats=handles.visual_stats;
        handles.brain=handles.visual_image;
        handles.t2star=handles.visual_t2star; %germuska@fmrib:2012 added to load appropriate t2star weigthed data
   
     case 2 
        handles.DataSet=1;
        handles.dzstats=handles.front_stats;
        handles.brain=handles.front_image;
        handles.t2star=handles.front_t2star;
end

%because images can differ in size reset voxel marker to 1,1
handles.x=1;
handles.y=1;

%update display with voxel location
stringX=num2str(handles.x); 
set(handles.x_text,'String',stringX);
stringY=num2str(handles.y); 
set(handles.y_text,'String',stringY);

% Update handles structure
guidata(hObject, handles);

%re-calculate z-stats
ScanButton_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function data_select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bandwidth_box_Callback(hObject, eventdata, handles)
% hObject    handle to bandwidth_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bandwidth_box as text
%        str2double(get(hObject,'String')) returns contents of bandwidth_box as a double

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
function bandwidth_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bandwidth_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in phase_menu.
function phase_menu_Callback(hObject, eventdata, handles)
% hObject    handle to phase_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns phase_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from phase_menu

val = get(hObject, 'Value');
switch val;
    case 1
        handles.phase='AP';
    case 2 
        handles.phase='LR';
end

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function phase_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phase_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in blip_menu.
function blip_menu_Callback(hObject, eventdata, handles)
% hObject    handle to blip_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns blip_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from blip_menu

val = get(hObject, 'Value');
switch val;
    case 1
        handles.blip=true;
    case 2 
        handles.blip=false;
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function blip_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blip_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in noise_box.
function noise_box_Callback(hObject, eventdata, handles)
% hObject    handle to noise_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noise_box

val=get(hObject,'Value'); 

switch val;
    case 0
        handles.regress=false;
    case 1 
        handles.regress=true;
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in im_weight_pop.
function im_weight_pop_Callback(hObject, eventdata, handles)
% hObject    handle to im_weight_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns im_weight_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from im_weight_pop

%if dataset exists update the display
if handles.DataReady>0
    plotOverlay(handles)
end


% --- Executes during object creation, after setting all properties.
function im_weight_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to im_weight_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
