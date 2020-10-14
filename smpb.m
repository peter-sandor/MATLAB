function varargout = smpb(varargin)
% SMPB M-file for smpb.fig
%      SMPB, by itself, creates a new SMPB or raises the existing
%      singleton*.
%
%      H = SMPB returns the handle to a new SMPB or the handle to
%      the existing singleton*.
%
%      SMPB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SMPB.M with the given input arguments.
%
%      SMPB('Property','Value',...) creates a new SMPB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before smpb_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to smpb_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help smpb

% Last Modified by GUIDE v2.5 03-Feb-2015 15:04:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @smpb_OpeningFcn, ...
                   'gui_OutputFcn',  @smpb_OutputFcn, ...
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


% --- Executes just before smpb is made visible.
function smpb_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to smpb (see VARARGIN)

% Choose default command line output for smpb
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes smpb wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = smpb_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load_trace.
function load_trace_Callback(hObject, eventdata, handles)
% hObject    handle to load_trace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile({'*.*','All Files (*.*)'},'Pick a trace');
temp=importdata(strcat(PathName,FileName));
if length(size(temp))==1
    handles.trace=temp;
elseif length(size(temp))==2
    handles.trace=squeeze(mean(temp,2));
end
guidata(hObject,handles)
plot(handles.axes1,handles.trace)



% --- Executes on button press in add_peak.
function add_peak_Callback(hObject, eventdata, handles)
% hObject    handle to add_peak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'peakdata')
    N=length(handles.peakdata);
else N=0;
end
handles.peakdata{N+1}.name=get(handles.Name,'String');
handles.peakdata{N+1}.index=[str2num(get(handles.peak_start,'String')) str2num(get(handles.peak_end,'String')) str2num(get(handles.bkg_start,'String')) str2num(get(handles.bkg_end,'String'))];
for ind1=1:N+1
    temp{ind1}=[handles.peakdata{ind1}.name ' (' num2str(handles.peakdata{ind1}.index) ')'];
end
set(handles.peaknames,'String',temp);
set(handles.peaknames,'Value',N+1);
guidata(hObject,handles)

% --- Executes on button press in zoom.
function zoom_Callback(hObject, eventdata, handles)
% hObject    handle to zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.popupmenu,'Value')==1
    zoom_start=max(1,str2num(get(handles.zoom_start,'String')));
    zoom_end=min(length(handles.trace),abs(str2num(get(handles.zoom_end,'String'))));
    plot(handles.axes1,zoom_start:zoom_end,handles.trace(zoom_start:zoom_end));
elseif get(handles.popupmenu,'Value')==2 && isfield(handles,'massaxis')
    mass_start=max(0,round(str2num(get(handles.zoom_start,'String'))));
    mass_end=round(min(handles.massaxis(end),str2num(get(handles.zoom_end,'String'))));
    temp=1:length(handles.massaxis);
    zoom_start=round(mean(temp(mass_start==round(handles.massaxis))));
    zoom_end=round(mean(temp(mass_end==round(handles.massaxis))));
    plot(handles.axes1,handles.massaxis(zoom_start:zoom_end),handles.trace(zoom_start:zoom_end));
end

function zoom_start_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zoom_start as text
%        str2double(get(hObject,'String')) returns contents of zoom_start as a double


% --- Executes during object creation, after setting all properties.
function zoom_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoom_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zoom_end_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zoom_end as text
%        str2double(get(hObject,'String')) returns contents of zoom_end as a double


% --- Executes during object creation, after setting all properties.
function zoom_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoom_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function name_Callback(hObject, eventdata, handles)
% hObject    handle to Name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Name as text
%        str2double(get(hObject,'String')) returns contents of Name as a double


% --- Executes during object creation, after setting all properties.
function name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function peak_start_Callback(hObject, eventdata, handles)
% hObject    handle to Peak_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Peak_start as text
%        str2double(get(hObject,'String')) returns contents of Peak_start as a double


% --- Executes during object creation, after setting all properties.
function peak_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Peak_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function peak_end_Callback(hObject, eventdata, handles)
% hObject    handle to peak_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of peak_end as text
%        str2double(get(hObject,'String')) returns contents of peak_end as a double


% --- Executes during object creation, after setting all properties.
function peak_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to peak_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function bkg_start_Callback(hObject, eventdata, handles)
% hObject    handle to bkg_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bkg_start as text
%        str2double(get(hObject,'String')) returns contents of bkg_start as a double


% --- Executes during object creation, after setting all properties.
function bkg_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bkg_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bkg_end_Callback(hObject, eventdata, handles)
% hObject    handle to bkg_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bkg_end as text
%        str2double(get(hObject,'String')) returns contents of bkg_end as a double


% --- Executes during object creation, after setting all properties.
function bkg_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bkg_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in peaknames.
function peaknames_Callback(hObject, eventdata, handles)
% hObject    handle to peaknames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns peaknames contents as cell array
%        contents{get(hObject,'Value')} returns selected item from peaknames


% --- Executes during object creation, after setting all properties.
function peaknames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to peaknames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in save_bins.
function save_bins_Callback(hObject, eventdata, handles)
% hObject    handle to save_bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uiputfile('TOF_bins.txt','Save file name');
fid=fopen([path file],'w');
for ind1=1:length(handles.peakdata)
    fprintf(fid,'%s \n',[handles.peakdata{ind1}.name ' ' num2str(handles.peakdata{ind1}.index)]);
end
fclose(fid);


% --- Executes on button press in load_bins.
function load_bins_Callback(hObject, eventdata, handles)
% hObject    handle to load_bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile({'*.*','All Files (*.*)'},'Pick a bin file.');
fid=fopen([PathName FileName]);
temp = textscan(fid, '%s %d %d %d %d');
for ind1=1:length(temp{1})
    handles.peakdata{ind1}.name=temp{1}{ind1};
    handles.peakdata{ind1}.index=[temp{2}(ind1) temp{3}(ind1) temp{4}(ind1) temp{5}(ind1)];
    temp2{ind1}=[handles.peakdata{ind1}.name ' (' num2str(handles.peakdata{ind1}.index) ')'];
end
set(handles.peaknames,'String',temp2);
guidata(hObject,handles);


% --- Executes on button press in delete_peak.
function delete_peak_Callback(hObject, eventdata, handles)
% hObject    handle to delete_peak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index_del=get(handles.peaknames,'Value');
handles.peakdata(index_del)=[];
N=length(handles.peakdata);
if N~=0
    for ind1=1:length(handles.peakdata)
        temp{ind1}=[handles.peakdata{ind1}.name ' (' num2str(handles.peakdata{ind1}.index) ')'];
    end
    set(handles.peaknames,'String',temp);
    set(handles.peaknames,'Value',N)
else set(handles.peaknames,'String','No peaks added');
    set(handles.peaknames,'Value',1)
end
guidata(hObject,handles);


% --- Executes on button press in calibrate.
function calibrate_Callback(hObject, eventdata, handles)
% hObject    handle to calibrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'trace')
    [handles.massaxis,handles.t0,handles.A]=TOFcalib(str2num(get(handles.index_mass1(2),'String')),str2num(get(handles.mass1(2),'String')),str2num(get(handles.index_mass2(2),'String')),str2num(get(handles.mass2(2),'String')),1:length(handles.trace));
    plot(handles.axes1,handles.massaxis,handles.trace)
else [handles.massaxis,handles.t0,handles.A]=TOFcalib(str2num(get(handles.index_mass1(2),'String')),str2num(get(handles.mass1(2),'String')),str2num(get(handles.index_mass2(2),'String')),str2num(get(handles.mass2(2),'String')),1:8000);
end
set(handles.cnst_t0,'String',num2str(handles.t0));
set(handles.cnst_A,'String',num2str(handles.A));
guidata(hObject,handles);


function index_mass1_Callback(hObject, eventdata, handles)
% hObject    handle to index_mass1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of index_mass1 as text
%        str2double(get(hObject,'String')) returns contents of index_mass1 as a double


% --- Executes during object creation, after setting all properties.
function index_mass1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to index_mass1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mass1_Callback(hObject, eventdata, handles)
% hObject    handle to mass1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mass1 as text
%        str2double(get(hObject,'String')) returns contents of mass1 as a double


% --- Executes during object creation, after setting all properties.
function mass1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mass1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function index_mass2_Callback(hObject, eventdata, handles)
% hObject    handle to index_mass2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of index_mass2 as text
%        str2double(get(hObject,'String')) returns contents of index_mass2 as a double


% --- Executes during object creation, after setting all properties.
function index_mass2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to index_mass2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mass2_Callback(hObject, eventdata, handles)
% hObject    handle to mass2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mass2 as text
%        str2double(get(hObject,'String')) returns contents of mass2 as a double


% --- Executes during object creation, after setting all properties.
function mass2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mass2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cnst_t0_Callback(hObject, eventdata, handles)
% hObject    handle to cnst_t0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cnst_t0 as text
%        str2double(get(hObject,'String')) returns contents of cnst_t0 as a double


% --- Executes during object creation, after setting all properties.
function cnst_t0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cnst_t0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cnst_A_Callback(hObject, eventdata, handles)
% hObject    handle to cnst_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cnst_A as text
%        str2double(get(hObject,'String')) returns contents of cnst_A as a double


% --- Executes during object creation, after setting all properties.
function cnst_A_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cnst_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in replot.
function replot_Callback(hObject, eventdata, handles)
% hObject    handle to replot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot(handles.axes1,handles.trace)



% --- Executes on selection change in popupmenu.
function popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu


% --- Executes during object creation, after setting all properties.
function popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in save_calib.
function save_calib_Callback(hObject, eventdata, handles)
% hObject    handle to save_calib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'massaxis')
    [file,path]=uiputfile('mass_axis.txt','Save mass axis');
    dlmwrite([path file],handles.massaxis,'\t');
    [file,path] = uiputfile('mass_calib_params.txt','Save calibration parameters');
    fid=fopen([path file],'w');
    fprintf(fid,'%s \n',['t0 ' get(handles.cnst_t0,'string')]);
    fprintf(fid,'%s \n',['A ' get(handles.cnst_A,'string')]);
    fclose(fid);
end