function varargout = pulse_analyzer(varargin)
% PULSE_ANALYZER MATLAB code for pulse_analyzer.fig
%      PULSE_ANALYZER, by itself, creates a new PULSE_ANALYZER or raises the existing
%      singleton*.
%
%      H = PULSE_ANALYZER returns the handle to a new PULSE_ANALYZER or the handle to
%      the existing singleton*.
%
%      PULSE_ANALYZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PULSE_ANALYZER.M with the given input arguments.
%
%      PULSE_ANALYZER('Property','Value',...) creates a new PULSE_ANALYZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pulse_analyzer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pulse_analyzer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pulse_analyzer

% Last Modified by GUIDE v2.5 21-Nov-2014 23:36:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pulse_analyzer_OpeningFcn, ...
                   'gui_OutputFcn',  @pulse_analyzer_OutputFcn, ...
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


% --- Executes just before pulse_analyzer is made visible.
function pulse_analyzer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pulse_analyzer (see VARARGIN)

% Choose default command line output for pulse_analyzer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pulse_analyzer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pulse_analyzer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load_spectrum.
function load_spectrum_Callback(hObject, eventdata, handles)
% hObject    handle to load_spectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sp_omega=SpectrumConvert;
[omega0 sigma FWHM]=calc_stats([sp_omega(:,1) sp_omega(:,2).^2]);
handles.omega0=omega0;
handles.FWHM=roundP(FWHM,4);
plot(handles.spectrum1,sp_omega(:,1),sp_omega(:,2).^2,'r');
setfigP(handles.spectrum1)
xlabel(handles.spectrum1,'\omega [rad/fs]')
ylabel(handles.spectrum1,'Intensity (normalized)')
xlim(handles.spectrum1,[min(sp_omega(:,1)) max(sp_omega(:,1))]);
ylim(handles.spectrum1,[0 1]);
handles.sp_omega=sp_omega;
% handles.omega0=sum(handles.sp_omega(:,1).*handles.sp_omega(:,2).^2)/sum(handles.sp_omega(:,2).^2);
set(handles.text_omega0,'String',handles.omega0);
set(handles.text_FWHM,'String',handles.FWHM);
set(handles.text_lambda0,'String',300/handles.omega0*2*pi);
set(handles.text_dlambda,'String',300/(handles.omega0/2/pi)^2*handles.FWHM/2/pi);

handles.td_TL=simulate_pulse_shaper([handles.sp_omega(:,1) handles.sp_omega(:,2)],'1');
handles.mask=get(handles.spectral_mask,'String');
if isempty(handles.mask)
    handles.mask='1';
end
[handles.td_data handles.sd_masked]=simulate_pulse_shaper([handles.sp_omega(:,1) handles.sp_omega(:,2)],handles.mask);
handles.phase_td=unwrap(angle(handles.td_data(:,2)));
intensity_td=abs(handles.td_data(:,2)).^2/max(abs(handles.td_data(:,2)).^2);
[a b tau]=calc_stats([handles.td_data(:,1) intensity_td]);
handles.tau=roundP(tau,1);
set(handles.text_tau,'String',handles.tau);
text_td{1}='\tau [fs]';
text_td{2}='Intensity (normalized)';
text_td{3}='Phase [pi]';
text_td{4}='';
TL_ip=interp1(handles.td_TL(:,1),handles.td_TL(:,2),handles.td_data(:,1));
plotyyP2(handles.td_data(:,1),[abs(TL_ip).^2/max(abs(TL_ip).^2) abs(handles.td_data(:,2)).^2/max(abs(handles.td_data(:,2)).^2)],handles.phase_td,text_td,handles.timedomain1)

field1=handles.td_data(:,2);
t0=sum(handles.td_data(:,1).*abs(handles.td_data(:,2)))/sum(abs(handles.td_data(:,2)));
timeax=handles.td_data(:,1)-t0;
Nt=length(timeax);
timeax_step=abs(timeax(2)-timeax(1));
freqax=map2colvec(FourierAxis(timeax));%+omega0;
tauax=map2colvec(-Nt*timeax_step:timeax_step:(Nt-1)*timeax_step);

field2=[zeros([Nt 1]); field1; zeros([Nt 1])];
signal_SHG=zeros([Nt 2*Nt]);
signal_SD=zeros([Nt 2*Nt]);

for ind1=1:2*Nt
    field2_shifted=circshift(field2,[-Nt+ind1 1]);
    signal_SHG(:,ind1)=map2rowvec(abs(fft(field1.^1.*field2_shifted(Nt+1:2*Nt).^1)).^2);
    signal_SD(:,ind1)=map2rowvec(abs(fft(field1.^2.*conj(field2_shifted(Nt+1:2*Nt)))).^2);
end

handles.signal_SHG=signal_SHG/max(max(signal_SHG));
handles.signal_SD=signal_SD/max(max(signal_SD));

imagescP(tauax,fftshift(freqax),fftshift(handles.signal_SHG,1),'Parent',handles.FROG_SHG)
xlabel(handles.FROG_SHG,'delay [fs]');
ylabel(handles.FROG_SHG,'\nu [PHz]');
title(handles.FROG_SHG,'SHG');
colormap bone;
CMAP=colormap;
CMAP=flipud(CMAP);
colormap(CMAP);

imagescP(tauax,fftshift(freqax),fftshift(handles.signal_SD,1),'Parent',handles.FROG_SD)
% xlabel(handles.FROG_SD,'t');
% ylabel(handles.FROG_SD,'\nu');
title(handles.FROG_SD,'SD');

% set(handles.text_tau,'String',0);
% plot(handles.timedomain1,[0],[0],'k');
% imagescP([0 1],[0 1],[0 0; 0 0],'Parent',handles.FROG_SHG)
% imagescP([0 1],[0 1],[0 0; 0 0],'Parent',handles.FROG_SD)
guidata(hObject,handles)



function spectral_mask_Callback(hObject, eventdata, handles)
% hObject    handle to spectral_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spectral_mask as text
%        str2double(get(hObject,'String')) returns contents of spectral_mask as a double


% --- Executes during object creation, after setting all properties.
function spectral_mask_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spectral_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calc_TD.
function calc_TD_Callback(hObject, eventdata, handles)
% hObject    handle to calc_TD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mask=get(handles.spectral_mask,'String');
if isempty(handles.mask)
    handles.mask='1';
end

[handles.td_data handles.sd_masked]=simulate_pulse_shaper([handles.sp_omega(:,1) handles.sp_omega(:,2)],handles.mask);
handles.phase_td=unwrap(angle(handles.td_data(:,2)));
intensity_td=abs(handles.td_data(:,2)).^2/max(abs(handles.td_data(:,2)).^2);
[a b tau]=calc_stats([handles.td_data(:,1) intensity_td]);
handles.tau=roundP(tau,1);
set(handles.text_tau,'String',handles.tau);
text_td{1}='\tau [fs]';
text_td{2}='Intensity (normalized)';
text_td{3}='Phase [pi]';
text_td{4}='';
TL_ip=interp1(handles.td_TL(:,1),handles.td_TL(:,2),handles.td_data(:,1));
plotyyP2(handles.td_data(:,1),[abs(TL_ip).^2/max(abs(TL_ip).^2) abs(handles.td_data(:,2)).^2/max(abs(handles.td_data(:,2)).^2)],handles.phase_td,text_td,handles.timedomain1)

omega=handles.sp_omega(:,1);
omega0=handles.omega0;
eval(['handles.phase_sp=unwrap(angle(' handles.mask '));']);
if length(handles.phase_sp)==1
    handles.phase_sp=handles.phase_sp*ones(size(handles.sp_omega(:,1)));
end
text_sp{1}='\omega [rad/fs]';
text_sp{2}='Intensity (normalized)';
text_sp{3}='Phase [pi]';
text_sp{4}='';
plotyyP2(handles.sp_omega(:,1),[abs(handles.sp_omega(:,2)).^2/max(abs(handles.sp_omega(:,2)).^2) abs(handles.sd_masked(:,2)).^2/max(abs(handles.sd_masked(:,2)).^2)],handles.phase_sp/pi,text_sp,handles.spectrum1)

field1=handles.td_data(:,2);
t0=sum(handles.td_data(:,1).*abs(handles.td_data(:,2)))/sum(abs(handles.td_data(:,2)));
timeax=handles.td_data(:,1)-t0;
Nt=length(timeax);
timeax_step=abs(timeax(2)-timeax(1));
freqax=map2colvec(FourierAxis(timeax));%+omega0;
tauax=map2colvec(-Nt*timeax_step:timeax_step:(Nt-1)*timeax_step);

field2=[zeros([Nt 1]); field1; zeros([Nt 1])];
signal_SHG=zeros([Nt 2*Nt]);
signal_SD=zeros([Nt 2*Nt]);

for ind1=1:2*Nt
    field2_shifted=circshift(field2,[-Nt+ind1 1]);
    signal_SHG(:,ind1)=map2rowvec(abs(fft(field1.^1.*field2_shifted(Nt+1:2*Nt).^1)).^2);
    signal_SD(:,ind1)=map2rowvec(abs(fft(field1.^2.*conj(field2_shifted(Nt+1:2*Nt)))).^2);
end

handles.signal_SHG=signal_SHG/max(max(signal_SHG));
handles.signal_SD=signal_SD/max(max(signal_SD));

imagescP(tauax,fftshift(freqax),fftshift(handles.signal_SHG,1),'Parent',handles.FROG_SHG)
xlabel(handles.FROG_SHG,'delay [fs]');
ylabel(handles.FROG_SHG,'\nu [PHz]');
title(handles.FROG_SHG,'SHG');
colormap bone;
CMAP=colormap;
CMAP=flipud(CMAP);
colormap(CMAP);

imagescP(tauax,fftshift(freqax),fftshift(handles.signal_SD,1),'Parent',handles.FROG_SD)
% xlabel(handles.FROG_SD,'t');
% ylabel(handles.FROG_SD,'\nu');
title(handles.FROG_SD,'SD');

guidata(hObject,handles);


function text_omega0_Callback(hObject, eventdata, handles)
% hObject    handle to text_omega0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_omega0 as text
%        str2double(get(hObject,'String')) returns contents of text_omega0 as a double


% --- Executes during object creation, after setting all properties.
function text_omega0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_omega0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_FWHM_Callback(hObject, eventdata, handles)
% hObject    handle to text_FWHM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_FWHM as text
%        str2double(get(hObject,'String')) returns contents of text_FWHM as a double


% --- Executes during object creation, after setting all properties.
function text_FWHM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_FWHM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_tau_Callback(hObject, eventdata, handles)
% hObject    handle to text_tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_tau as text
%        str2double(get(hObject,'String')) returns contents of text_tau as a double


% --- Executes during object creation, after setting all properties.
function text_tau_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_lambda0_Callback(hObject, eventdata, handles)
% hObject    handle to text_lambda0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_lambda0 as text
%        str2double(get(hObject,'String')) returns contents of text_lambda0 as a double


% --- Executes during object creation, after setting all properties.
function text_lambda0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_lambda0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_dlambda_Callback(hObject, eventdata, handles)
% hObject    handle to text_dlambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_dlambda as text
%        str2double(get(hObject,'String')) returns contents of text_dlambda as a double


% --- Executes during object creation, after setting all properties.
function text_dlambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_dlambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
