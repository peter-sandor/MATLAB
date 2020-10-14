function varargout = popup_lns(varargin)
% POPUP_LNS MATLAB code for popup_lns.fig
%      POPUP_LNS, by itself, creates a new POPUP_LNS or raises the existing
%      singleton*.
%
%      H = POPUP_LNS returns the handle to a new POPUP_LNS or the handle to
%      the existing singleton*.
%
%      POPUP_LNS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POPUP_LNS.M with the given input arguments.
%
%      POPUP_LNS('Property','Value',...) creates a new POPUP_LNS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before popup_lns_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to popup_lns_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help popup_lns

% Last Modified by GUIDE v2.5 06-Jan-2014 22:06:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @popup_lns_OpeningFcn, ...
                   'gui_OutputFcn',  @popup_lns_OutputFcn, ...
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


% --- Executes just before popup_lns is made visible.
function popup_lns_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to popup_lns (see VARARGIN)

% Choose default command line output for popup_lns
handles.output = [];
handles.fighandle=hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes popup_lns wait for user response (see UIRESUME)
uiwait(handles.fighandle);


% --- Outputs from this function are returned to the command line.
function varargout = popup_lns_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
close(handles.fighandle);


function fsp_param_Callback(hObject, eventdata, handles)
% hObject    handle to fsp_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fsp_param as text
%        str2double(get(hObject,'String')) returns contents of fsp_param as a double
handles.lastdata=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function fsp_param_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fsp_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Done_button.
function Done_button_Callback(hObject, eventdata, handles)
% hObject    handle to Done_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output=handles.lastdata;
guidata(hObject, handles);
uiresume(handles.fighandle);
