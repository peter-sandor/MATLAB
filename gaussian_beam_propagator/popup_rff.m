function varargout = popup_rff(varargin)
% POPUP_RFF MATLAB code for popup_rff.fig
%      POPUP_RFF, by itself, creates a new POPUP_RFF or raises the existing
%      singleton*.
%
%      H = POPUP_RFF returns the handle to a new POPUP_RFF or the handle to
%      the existing singleton*.
%
%      POPUP_RFF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POPUP_RFF.M with the given input arguments.
%
%      POPUP_RFF('Property','Value',...) creates a new POPUP_RFF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before popup_rff_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to popup_rff_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help popup_rff

% Last Modified by GUIDE v2.5 06-Jan-2014 22:13:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @popup_rff_OpeningFcn, ...
                   'gui_OutputFcn',  @popup_rff_OutputFcn, ...
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


% --- Executes just before popup_rff is made visible.
function popup_rff_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to popup_rff (see VARARGIN)

% Choose default command line output for popup_rff
handles.output = [];
handles.fighandle=hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes popup_rff wait for user response (see UIRESUME)
uiwait(handles.fighandle);


% --- Outputs from this function are returned to the command line.
function varargout = popup_rff_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
close(handles.fighandle);


function n1_param_Callback(hObject, eventdata, handles)
% hObject    handle to n1_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n1_param as text
%        str2double(get(hObject,'String')) returns contents of n1_param as a double
handles.n1=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function n1_param_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n1_param (see GCBO)
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
handles.output=[handles.n1 handles.n2];
guidata(hObject, handles);
uiresume(handles.fighandle);



function n2_param_Callback(hObject, eventdata, handles)
% hObject    handle to n2_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n2_param as text
%        str2double(get(hObject,'String')) returns contents of n2_param as a double
handles.n2=str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function n2_param_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n2_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
