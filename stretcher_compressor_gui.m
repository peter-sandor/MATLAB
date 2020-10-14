function varargout = stretcher_compressor_gui(varargin)
% STRETCHER_COMPRESSOR_GUI MATLAB code for stretcher_compressor_gui.fig
%      STRETCHER_COMPRESSOR_GUI, by itself, creates a new STRETCHER_COMPRESSOR_GUI or raises the existing
%      singleton*.
%
%      H = STRETCHER_COMPRESSOR_GUI returns the handle to a new STRETCHER_COMPRESSOR_GUI or the handle to
%      the existing singleton*.
%
%      STRETCHER_COMPRESSOR_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STRETCHER_COMPRESSOR_GUI.M with the given input arguments.
%
%      STRETCHER_COMPRESSOR_GUI('Property','Value',...) creates a new STRETCHER_COMPRESSOR_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stretcher_compressor_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stretcher_compressor_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stretcher_compressor_gui

% Last Modified by GUIDE v2.5 21-May-2014 18:52:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stretcher_compressor_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @stretcher_compressor_gui_OutputFcn, ...
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


% --- Executes just before stretcher_compressor_gui is made visible.
function stretcher_compressor_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stretcher_compressor_gui (see VARARGIN)

% Choose default command line output for stretcher_compressor_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes stretcher_compressor_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = stretcher_compressor_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function lambda_Callback(hObject, eventdata, handles)
% hObject    handle to lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda as text
%        str2double(get(hObject,'String')) returns contents of lambda as a double


% --- Executes during object creation, after setting all properties.
function lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gc_Callback(hObject, eventdata, handles)
% hObject    handle to Gc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gc as text
%        str2double(get(hObject,'String')) returns contents of Gc as a double


% --- Executes during object creation, after setting all properties.
function Gc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gamma_c_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gamma_c as text
%        str2double(get(hObject,'String')) returns contents of gamma_c as a double


% --- Executes during object creation, after setting all properties.
function gamma_c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gs_Callback(hObject, eventdata, handles)
% hObject    handle to Gs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gs as text
%        str2double(get(hObject,'String')) returns contents of Gs as a double


% --- Executes during object creation, after setting all properties.
function Gs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gamma_s_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gamma_s as text
%        str2double(get(hObject,'String')) returns contents of gamma_s as a double


% --- Executes during object creation, after setting all properties.
function gamma_s_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
temp=get(handles.gamma_s,'String');
if iscell(temp)
    eval(['gamma_s=' temp{1} ';']);
else
    eval(['gamma_s=' temp ';']);
end
% handles.gamma_s=gamma_s;

temp=get(handles.gamma_c,'String');
if iscell(temp)
    eval(['gamma_c=' temp{1} ';']);
else
    eval(['gamma_c=' temp ';']);
end
% handles.gamma_c=gamma_c;

temp=get(handles.Gc,'String');
if iscell(temp)
    eval(['Gc=' temp{1} ';']);
else
    eval(['Gc=' temp ';']);
end
% handles.Gc=Gc;

temp=get(handles.Gs,'String');
if iscell(temp)
    eval(['Gs=' temp{1} ';']);
else
    eval(['Gs=' temp ';']);
end
% handles.Gs=Gs;

temp=get(handles.lambda,'String');
eval(['lambda=' temp ';']);
% handles.lambda=lambda;

% [phi2,phi3]=calc_stretcher_compressor(handles.Gc,handles.gamma_c,handles.Gs,handles.gamma_s,handles.lambda);
% handles.N=[length(handles.Gc) length(handles.gamma_c) length(handles.Gs) length(handles.gamma_s)];
[phi2,phi3]=calc_stretcher_compressor(Gc,gamma_c,Gs,gamma_s,lambda);
temp=get(handles.mat2,'String');
mat2=str2num(temp);
temp=get(handles.mat3,'String');
mat3=str2num(temp);
handles.N=[length(Gc) length(gamma_c) length(Gs) length(gamma_s)];
N3=length(handles.N(handles.N>1));

if N3==1
    index1=(1:4).*(handles.N>1);
    index1(index1==0)=[];
    switch index1
        case 1
            x=Gc;
            text_xlabel='G_c';
        case 2
            x=gamma_c;
            text_xlabel='\gamma_c';
        case 3
            x=Gs;
            text_xlabel='l_{eff}';
        case 4
            x=gamma_s;
            text_xlabel='\gamma_s';
    end
    plot(handles.phi2,x,phi2,'k')
    setfigP(handles.phi2)
    xlabel(handles.phi2,text_xlabel);
    ylabel(handles.phi2,'\Phi_2 [fs^2]');
    plot(handles.phi3,x,phi3,'k')
    setfigP(handles.phi3)
	xlabel(handles.phi3,text_xlabel);
    ylabel(handles.phi3,'\Phi_3 [fs^3]');
elseif N3==2
    N2=num2str(handles.N>1);
    N2(ismember(N2,' '))=[];
    switch num2str(map2rowvec(N2),'1.0d')
        case '1100'
            x=gamma_c;
            y=Gc;
            text_xlabel='\gamma_c';
            text_ylabel='Gc';
        case '1010'    
            x=Gs;
            y=Gc;
            text_xlabel='l_{eff}';
            text_ylabel='Gc';
        case '1001'
            x=gamma_s;
            y=Gc;
            text_xlabel='\gamma_s';
            text_ylabel='Gc';
        case '0110'
            x=Gs;
            y=gamma_c;
            text_xlabel='l_{eff}';
            text_ylabel='\gamma_c';
        case '0101'
            x=gamma_s;
            y=gamma_c;
            text_xlabel='\gamma_s';
            text_ylabel='\gamma_c';
        case '0011'
            x=gamma_s;
            y=Gs;
            text_xlabel='\gamma_s';
            text_ylabel='l_{eff}';
    end

    if get(handles.select_plot,'Value')==1
        colormap jet;
        imagescP(x,y,phi2+mat2,'Parent',handles.phi2);
        xlabel(handles.phi2,text_xlabel);
        ylabel(handles.phi2,text_ylabel);
        title(handles.phi2,'2nd order')
        colorbar('peer',handles.phi2,'Eastoutside')
        imagescP(x,y,phi3+mat3,'Parent',handles.phi3);
        xlabel(handles.phi3,text_xlabel);
        ylabel(handles.phi3,text_ylabel);
        title(handles.phi3,'3rd order')
        colorbar('peer',handles.phi3,'Eastoutside')
    elseif get(handles.select_plot,'Value')==2
        temp=get(handles.thrs2,'String');
        thrs2=str2num(temp);
        temp=get(handles.thrs3,'String');
        thrs3=str2num(temp);
        img_th2=(abs(phi2+mat2)<=thrs2);
        img_th3=(abs(phi3+mat3)<=thrs3);
        colormap bone;
        colormap(flipud(colormap));
        imagescP(x,y,img_th2+0.2*img_th3,'Parent',handles.phi2);
        xlabel(handles.phi2,text_xlabel);
        ylabel(handles.phi2,text_ylabel);
        title(handles.phi2,'2nd order')
        imagescP(x,y,img_th3+0.2*img_th2,'Parent',handles.phi3);
        xlabel(handles.phi3,text_xlabel);
        ylabel(handles.phi3,text_ylabel);
        title(handles.phi3,'3rd order')
    end
	set(handles.phi2,'XGrid','on','YGrid','on');
    set(handles.phi3,'XGrid','on','YGrid','on');
end
guidata(hObject,handles);


% --- Executes on selection change in select_plot.
function select_plot_Callback(hObject, eventdata, handles)
% hObject    handle to select_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns select_plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from select_plot


% --- Executes during object creation, after setting all properties.
function select_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thrs2_Callback(hObject, eventdata, handles)
% hObject    handle to thrs2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thrs2 as text
%        str2double(get(hObject,'String')) returns contents of thrs2 as a double


% --- Executes during object creation, after setting all properties.
function thrs2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thrs2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thrs3_Callback(hObject, eventdata, handles)
% hObject    handle to thrs3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thrs3 as text
%        str2double(get(hObject,'String')) returns contents of thrs3 as a double


% --- Executes during object creation, after setting all properties.
function thrs3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thrs3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mat2_Callback(hObject, eventdata, handles)
% hObject    handle to mat2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mat2 as text
%        str2double(get(hObject,'String')) returns contents of mat2 as a double


% --- Executes during object creation, after setting all properties.
function mat2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mat2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mat3_Callback(hObject, eventdata, handles)
% hObject    handle to mat3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mat3 as text
%        str2double(get(hObject,'String')) returns contents of mat3 as a double


% --- Executes during object creation, after setting all properties.
function mat3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mat3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
