function varargout = gaussbeam(varargin)
% GAUSSBEAM M-file for gaussbeam.fig
%      GAUSSBEAM, by itself, creates a new GAUSSBEAM or raises the existing
%      singleton*.
%
%      H = GAUSSBEAM returns the handle to a new GAUSSBEAM or the handle to
%      the existing singleton*.
%
%      GAUSSBEAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GAUSSBEAM.M with the given input arguments.
%
%      GAUSSBEAM('Property','Value',...) creates a new GAUSSBEAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before smpb_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gaussbeam_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gaussbeam

% Last Modified by GUIDE v2.5 16-Apr-2018 11:07:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gaussbeam_OpeningFcn, ...
                   'gui_OutputFcn',  @gaussbeam_OutputFcn, ...
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


% --- Executes just before gaussbeam is made visible.
function gaussbeam_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gaussbeam (see VARARGIN)

% Choose default command line output for gaussbeam
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gaussbeam wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gaussbeam_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Calculate.
function Calculate_Callback(hObject, eventdata, handles)
% hObject    handle to Calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles.M2=12;
temp=get(handles.input_M2,'String');
M=sqrt(str2double(temp(2)));
handles.input.M=M;
temp=get(handles.input_w0,'String');
handles.input.w=str2double(temp(2))/M; % [mm]
temp=get(handles.input_R0,'String');
handles.input.R=str2double(temp(2)); % [mm]
temp=get(handles.input_lambda,'String');
handles.input.lambda=str2double(temp(2))*1e-6; % [mm]

Nplot=200;
handles.ind_fsp=[];
handles.ind_lns=[];
handles.ind_rff=[];
handles.ind_rfc=[];
for ind1=1:length(handles.system)
    if strcmp(handles.system(ind1).element,'fsp')
        handles.ind_fsp=cat(2,handles.ind_fsp,ind1);
    elseif strcmp(handles.system(ind1).element,'lns')
        handles.ind_lns=cat(2,handles.ind_lns,ind1);
    elseif strcmp(handles.system(ind1).element,'rff')
        handles.ind_rff=cat(2,handles.ind_rff,ind1);
    elseif strcmp(handles.system(ind1).element,'rfc')
        handles.ind_rfc=cat(2,handles.ind_rfc,ind1);  
    end
end
handles.zplot=[];
handles.w=[];
handles.R=[];
ind3=1;

for ind2=handles.ind_fsp
%     z_start=sum(system.s(ind_fsp(ind_fsp<ind2)));
    z_start=0;
    z_end=handles.system(ind2).param;
    z=z_start:(z_end-z_start)/(Nplot-1):z_end;
    if isempty(handles.ind_fsp(handles.ind_fsp<ind2))
        handles.zplot=cat(2,handles.zplot,z);
    else
        handles.zplot=cat(2,handles.zplot,sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind2)).param])+z);
    end
    if sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind2)).param])<str2double(get(handles.z_value,'String')) && sum([handles.system(handles.ind_fsp(handles.ind_fsp<=ind2)).param])>str2double(get(handles.z_value,'String'))
        subsystem2=handles.system(1:ind2);
        subsystem2(ind2).param=str2double(get(handles.z_value,'String'))-sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind2)).param]);
        output2=beam_tracer(handles.input,subsystem2);
        handles.calced=output2;
        waist_loc = str2double(get(handles.z_value,'String')) - handles.calced.z_waist;
        set(handles.w_out,'string',num2str(handles.calced.w*M));
        set(handles.R_out,'string',num2str(roundP(handles.calced.R,1)));
        set(handles.lambda_out,'string',num2str(handles.calced.lambda*1e6));
        set(handles.w0,'string',num2str(handles.calced.beam_waist*M));
        set(handles.rayleigh,'string',num2str(roundP(handles.calced.beam_waist^2*pi/handles.calced.lambda,1)));
        set(handles.z_waist,'string',num2str(roundP(waist_loc,1)));
    end
%     input.element=system.element(1:ind2,:);
    subsystem=handles.system(1:ind2);
    for ind1=1:Nplot
        subsystem(ind2).param=z(ind1);
        output=beam_tracer(handles.input,subsystem);
        handles.R(ind3)=output.R;
        handles.w(ind3)=output.w*M;
        ind3=ind3+1;
    end
end
guidata(hObject,handles)

if get(handles.popupmenu,'value')==1 % plot beam size
    plot(handles.axes1,handles.zplot,handles.w,'k')
    xlabel('z [mm]')
    ylabel('M\cdotw [mm]')
    zrange=max(handles.zplot)-min(handles.zplot);
    wmax=1.3*max(handles.w);
    for ind1=handles.ind_lns
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
        temp=1.3*handles.w(handles.zplot==dist);
        line([dist dist],[0 wmax],'color','r');
        if handles.system(ind1).param>0
            line([dist-0.03*zrange dist],[0.95*wmax wmax],'color','r');
            line([dist dist+0.03*zrange],[wmax 0.95*wmax],'color','r');
        elseif handles.system(ind1).param<0
            line([dist-0.03*zrange dist],[1.05*wmax wmax],'color','r');
            line([dist dist+0.03*zrange],[wmax 1.05*wmax],'color','r');
        end
    end

    for ind1=handles.ind_rff
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
    %     temp=1.3*handles.w(handles.zplot==dist);
        line([dist dist],[0 wmax],'color','b');
        if handles.system(ind1).param(2)/handles.system(ind1).param(1)<=1
            line([dist-0.03*zrange dist],[wmax wmax],'color','b');
        else
            line([dist dist+0.03*zrange],[wmax wmax],'color','b');
        end
    end

    for ind1=handles.ind_rfc
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
        temp=1.3*handles.w(handles.zplot==dist);
        line([dist dist],[0 wmax],'color','g');
        if handles.system(ind1).param(2)/handles.system(ind1).param(1)<=1
            if handles.system(ind1).param(3)>0
            line([dist-0.03*zrange dist],[1.05*wmax wmax],'color','g');
            elseif handles.system(ind1).param(3)<=0
            line([dist-0.03*zrange dist],[0.95*wmax wmax],'color','g');
            end
        else
            if handles.system(ind1).param(3)>0
            line([dist dist+0.03*zrange],[wmax 0.95*wmax],'color','g');
            elseif handles.system(ind1).param(3)<=0
            line([dist dist+0.03*zrange],[wmax 1.05*wmax],'color','g');
            end
        end
    end
    
elseif get(handles.popupmenu,'value')==2 %plot curvature
    plot(handles.axes1,handles.zplot,handles.R,'k')
    xlabel('z [mm]')
    ylabel('R [mm]')
    zrange=max(handles.zplot)-min(handles.zplot);
    Rmax=1.3*max(handles.R);
    for ind1=handles.ind_lns
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
        temp=1.3*handles.R(handles.zplot==dist);
        line([dist dist],[0 Rmax],'color','r');
        if handles.system(ind1).param>0
            line([dist-0.03*zrange dist],[0.95*Rmax Rmax],'color','r');
            line([dist dist+0.03*zrange],[Rmax 0.95*Rmax],'color','r');
        elseif handles.system(ind1).param<0
            line([dist-0.03*zrange dist],[1.05*Rmax Rmax],'color','r');
            line([dist dist+0.03*zrange],[Rmax 1.05*Rmax],'color','r');
        end
    end

    for ind1=handles.ind_rff
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
    %     temp=1.3*handles.w(handles.zplot==dist);
        line([dist dist],[0 Rmax],'color','b');
        if handles.system(ind1).param(2)/handles.system(ind1).param(1)<=1
            line([dist-0.03*zrange dist],[Rmax Rmax],'color','b');
        else
            line([dist dist+0.03*zrange],[Rmax Rmax],'color','b');
        end
    end

    for ind1=handles.ind_rfc
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
        temp=1.3*handles.R(handles.zplot==dist);
        line([dist dist],[0 Rmax],'color','g');
        if handles.system(ind1).param(2)/handles.system(ind1).param(1)<=1
            if handles.system(ind1).param(3)>0
            line([dist-0.03*zrange dist],[1.05*Rmax Rmax],'color','g');
            elseif handles.system(ind1).param(3)<=0
            line([dist-0.03*zrange dist],[0.95*Rmax Rmax],'color','g');
            end
        else
            if handles.system(ind1).param(3)>0
            line([dist dist+0.03*zrange],[Rmax 0.95*Rmax],'color','g');
            elseif handles.system(ind1).param(3)<=0
            line([dist dist+0.03*zrange],[Rmax 1.05*Rmax],'color','g');
            end
        end
    end
end

% --- Executes on button press in plot_sep.
function plot_sep_Callback(hObject, eventdata, handles)
% hObject    handle to plot_sep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hfig1=figure;hax1=axes;
str_to_print{1}=['M = ' num2str(handles.input.M)];
str_to_print{2}=['M*w = ' num2str(handles.input.M*handles.input.w) ' mm'];
str_to_print{3}=['R = ' num2str(handles.input.R) ' mm'];
str_to_print{4}=['\lambda = ' num2str(handles.input.lambda*1e6) ' nm'];
if get(handles.popupmenu,'value')==1 % plot beam size
    plot(hax1,handles.zplot,handles.w,'k')
    xlabel('z [mm]')
    ylabel('M\cdotw [mm]')
    zrange=max(handles.zplot)-min(handles.zplot);
    wmax=1.3*max(handles.w);
    for ind1=handles.ind_lns
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
        temp=1.3*handles.w(handles.zplot==dist);
        line([dist dist],[0 wmax],'color','r');
        text(dist-0.01*zrange,1.1*wmax,num2str(handles.system(ind1).param))
        if handles.system(ind1).param>0
            line([dist-0.03*zrange dist],[0.95*wmax wmax],'color','r');
            line([dist dist+0.03*zrange],[wmax 0.95*wmax],'color','r');
        elseif handles.system(ind1).param<0
            line([dist-0.03*zrange dist],[1.05*wmax wmax],'color','r');
            line([dist dist+0.03*zrange],[wmax 1.05*wmax],'color','r');
        end
    end

    for ind1=handles.ind_rff
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
    %     temp=1.3*handles.w(handles.zplot==dist);
        line([dist dist],[0 wmax],'color','b');
        if handles.system(ind1).param(2)/handles.system(ind1).param(1)<=1
            line([dist-0.03*zrange dist],[wmax wmax],'color','b');
        else
            line([dist dist+0.03*zrange],[wmax wmax],'color','b');
        end
    end

    for ind1=handles.ind_rfc
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
        temp=1.3*handles.w(handles.zplot==dist);
        line([dist dist],[0 wmax],'color','g');
        if handles.system(ind1).param(2)/handles.system(ind1).param(1)<=1
            if handles.system(ind1).param(3)>0
            line([dist-0.03*zrange dist],[1.05*wmax wmax],'color','g');
            elseif handles.system(ind1).param(3)<=0
            line([dist-0.03*zrange dist],[0.95*wmax wmax],'color','g');
            end
        else
            if handles.system(ind1).param(3)>0
            line([dist dist+0.03*zrange],[wmax 0.95*wmax],'color','g');
            elseif handles.system(ind1).param(3)<=0
            line([dist dist+0.03*zrange],[wmax 1.05*wmax],'color','g');
            end
        end
    end
    
elseif get(handles.popupmenu,'value')==2 %plot curvature
    plot(hax1,handles.zplot,handles.R,'k')
    xlabel('z [mm]')
    ylabel('R [mm]')
    zrange=max(handles.zplot)-min(handles.zplot);
    Rmax=1.3*max(handles.R);
    for ind1=handles.ind_lns
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
        temp=1.3*handles.R(handles.zplot==dist);
        line([dist dist],[0 Rmax],'color','r');
        text(dist,1.1*Rmax,num2str(handles.system(ind1).param))
        if handles.system(ind1).param>0
            line([dist-0.03*zrange dist],[0.95*Rmax Rmax],'color','r');
            line([dist dist+0.03*zrange],[Rmax 0.95*Rmax],'color','r');
        elseif handles.system(ind1).param<0
            line([dist-0.03*zrange dist],[1.05*Rmax Rmax],'color','r');
            line([dist dist+0.03*zrange],[Rmax 1.05*Rmax],'color','r');
        end
    end

    for ind1=handles.ind_rff
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
    %     temp=1.3*handles.w(handles.zplot==dist);
        line([dist dist],[0 Rmax],'color','b');
        if handles.system(ind1).param(2)/handles.system(ind1).param(1)<=1
            line([dist-0.03*zrange dist],[Rmax Rmax],'color','b');
        else
            line([dist dist+0.03*zrange],[Rmax Rmax],'color','b');
        end
    end

    for ind1=handles.ind_rfc
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
        temp=1.3*handles.R(handles.zplot==dist);
        line([dist dist],[0 Rmax],'color','g');
        if handles.system(ind1).param(2)/handles.system(ind1).param(1)<=1
            if handles.system(ind1).param(3)>0
            line([dist-0.03*zrange dist],[1.05*Rmax Rmax],'color','g');
            elseif handles.system(ind1).param(3)<=0
            line([dist-0.03*zrange dist],[0.95*Rmax Rmax],'color','g');
            end
        else
            if handles.system(ind1).param(3)>0
            line([dist dist+0.03*zrange],[Rmax 0.95*Rmax],'color','g');
            elseif handles.system(ind1).param(3)<=0
            line([dist dist+0.03*zrange],[Rmax 1.05*Rmax],'color','g');
            end
        end
    end
end
title(str_to_print);

% --- Executes on button press in add.
function add_Callback(hObject, eventdata, handles)
% hObject    handle to add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'system')
    N=length(handles.system);
else N=0;
end

if get(handles.typeselect,'Value')==1 % add freespace
    handles.system(N+1).element='fsp';
    handles.system(N+1).param=popup_fsp;
elseif get(handles.typeselect,'Value')==2 % add thin lens
    handles.system(N+1).element='lns';
    handles.system(N+1).param=popup_lns;
elseif get(handles.typeselect,'Value')==3 % add planar boundary
    handles.system(N+1).element='rff';
    handles.system(N+1).param=popup_rff;   
elseif get(handles.typeselect,'Value')==4 % add spherical boundary
	handles.system(N+1).element='rfc';
    handles.system(N+1).param=popup_rfc;
end

% handles.system(N+1).element=get(handles.Name,'String');
% handles.system(N+1).data=[str2num(get(handles.peak_start,'String')) str2num(get(handles.peak_end,'String')) str2num(get(handles.bkg_start,'String')) str2num(get(handles.bkg_end,'String'))];

for ind1=1:N+1
    temp{ind1}=[handles.system(ind1).element ' (' num2str(handles.system(ind1).param) ')'];
end
set(handles.elementlist,'String',temp);
set(handles.elementlist,'Value',N+1);
guidata(hObject,handles)

% --- Executes on button press in zoom.
function zoom_Callback(hObject, eventdata, handles)
% hObject    handle to zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zoom_start=max(0,str2num(get(handles.zoom_start,'String')));
zoom_end=min(max(handles.zplot),abs(str2num(get(handles.zoom_end,'String'))));
index_start=min(find(handles.zplot>zoom_start));
index_end=max(find(handles.zplot<zoom_end));
zmin=min(handles.zplot(index_start:index_end));
zmax=max(handles.zplot(index_start:index_end));
zrange=zmax-zmin;
if get(handles.popupmenu,'Value')==1
    plot(handles.axes1,handles.zplot(index_start:index_end),handles.w(index_start:index_end),'k');
    wmax=1.3*max(handles.w(index_start:index_end));
    for ind1=handles.ind_lns
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
        if dist>zmin && dist<zmax
            temp=1.3*handles.w(handles.zplot==dist);
            line([dist dist],[0 wmax],'color','r');
            if handles.system(ind1).param>0
                line([dist-0.03*zrange dist],[0.95*wmax wmax],'color','r');
                line([dist dist+0.03*zrange],[wmax 0.95*wmax],'color','r');
            elseif handles.system(ind1).param<0
                line([dist-0.03*zrange dist],[1.05*wmax wmax],'color','r');
                line([dist dist+0.03*zrange],[wmax 1.05*wmax],'color','r');
            end
        end
    end
    
    for ind1=handles.ind_rff
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
        if dist>zmin && dist<zmax
            temp=1.3*handles.w(handles.zplot==dist);
            line([dist dist],[0 wmax],'color','b');
            if handles.system(ind1).param(2)/handles.system(ind1).param(1)<=1
                line([dist-0.03*zrange dist],[wmax wmax],'color','b');
            else
                line([dist dist+0.03*zrange],[wmax wmax],'color','b');
            end
        end
    end

    for ind1=handles.ind_rfc
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
%         temp=1.3*handles.w(handles.zplot==dist);
        if dist>zmin && dist<zmax
            line([dist dist],[0 wmax],'color','g');
            if handles.system(ind1).param(2)/handles.system(ind1).param(1)<=1
                if handles.system(ind1).param(3)>0
                line([dist-0.03*zrange dist],[1.05*wmax wmax],'color','g');
                elseif handles.system(ind1).param(3)<=0
                line([dist-0.03*zrange dist],[0.95*wmax wmax],'color','g');
                end
            else
                if handles.system(ind1).param(3)>0
                line([dist dist+0.03*zrange],[wmax 0.95*wmax],'color','g');
                elseif handles.system(ind1).param(3)<=0
                line([dist dist+0.03*zrange],[wmax 1.05*wmax],'color','g');
                end
            end
        end
    end
elseif get(handles.popupmenu,'Value')==2
    plot(handles.axes1,handles.zplot(index_start:index_end),handles.R(index_start:index_end),'k');
    Rmax=1.3*max(handles.R(index_start:index_end));
    for ind1=handles.ind_lns
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
        if dist>zmin && dist<zmax
            temp=1.3*handles.w(handles.zplot==dist);
            line([dist dist],[0 Rmax],'color','r');
            if handles.system(ind1).param>0
                line([dist-0.03*zrange dist],[0.95*Rmax Rmax],'color','r');
                line([dist dist+0.03*zrange],[Rmax 0.95*Rmax],'color','r');
            elseif handles.system(ind1).param<0
                line([dist-0.03*zrange dist],[1.05*Rmax Rmax],'color','r');
                line([dist dist+0.03*zrange],[Rmax 1.05*Rmax],'color','r');
            end
        end
    end
    
    for ind1=handles.ind_rff
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
        if dist>zmin && dist<zmax
            temp=1.3*handles.R(handles.zplot==dist);
            line([dist dist],[0 Rmax],'color','b');
            if handles.system(ind1).param(2)/handles.system(ind1).param(1)<=1
                line([dist-0.03*zrange dist],[Rmax Rmax],'color','b');
            else
                line([dist dist+0.03*zrange],[Rmax Rmax],'color','b');
            end
        end
    end

    for ind1=handles.ind_rfc
        dist=sum([handles.system(handles.ind_fsp(handles.ind_fsp<ind1)).param]);
%         temp=1.3*handles.w(handles.zplot==dist);
        if dist>zmin && dist<zmax
            line([dist dist],[0 Rmax],'color','g');
            if handles.system(ind1).param(2)/handles.system(ind1).param(1)<=1
                if handles.system(ind1).param(3)>0
                line([dist-0.03*zrange dist],[1.05*Rmax Rmax],'color','g');
                elseif handles.system(ind1).param(3)<=0
                line([dist-0.03*zrange dist],[0.95*Rmax Rmax],'color','g');
                end
            else
                if handles.system(ind1).param(3)>0
                line([dist dist+0.03*zrange],[Rmax 0.95*Rmax],'color','g');
                elseif handles.system(ind1).param(3)<=0
                line([dist dist+0.03*zrange],[Rmax 1.05*Rmax],'color','g');
                end
            end
        end
    end
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




% --- Executes on selection change in elementlist.
function elementlist_Callback(hObject, eventdata, handles)
% hObject    handle to elementlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns elementlist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from elementlist


% --- Executes during object creation, after setting all properties.
function elementlist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elementlist (see GCBO)
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
[file,path] = uiputfile('system.txt','Save file name');
fid=fopen([path file],'w');
for ind1=1:length(handles.system)
    fprintf(fid,'%s \n',[handles.system(ind1).element ' ' num2str(handles.system(ind1).param)]);
end
fclose(fid);


% --- Executes on button press in load_bins.
function load_bins_Callback(hObject, eventdata, handles)
% hObject    handle to load_bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile({'*.*','All Files (*.*)'},'Pick a file.');
fid=fopen([PathName FileName]);
frewind(fid);
handles.system=[];
ind1=1;
while ~feof(fid)
    temp=fgetl(fid);
    if temp~=-1
        handles.system(ind1).element=temp(1:3);
        handles.system(ind1).param=str2num(temp(4:end));
    end
    ind1=ind1+1;
end
% % temp = textscan(fid, '%s %d %d %d');
% temp2=fgetl(fid);
% for ind1=1:length(temp{1})
%     handles.system(ind1).element=temp{1}{ind1};
%     if strcmp(handles.system(ind1).element,'fsp') % modify freespace
%         handles.system(ind1).param=[temp{2}(ind1)];;
%     elseif strcmp(handles.system(ind1).element,'lns') % modify thin lens
%         handles.system(ind1).param=[temp{2}(ind1)];
%     elseif strcmp(handles.system(ind1).element,'rff') % modify planar boundary
%         handles.system(ind1).param=[temp{2}(ind1) temp{3}(ind1)];  
%     elseif strcmp(handles.system(ind1).element,'rfc') % modify curved boundary
%         handles.system(ind1).param=[temp{2}(ind1) temp{3}(ind1) temp{4}(ind1)];
%     end
%     handles.system(ind1).param=[temp{2}(ind1) temp{3}(ind1) temp{4}(ind1)];
%     temp2{ind1}=[handles.system(ind1).element ' (' num2str(handles.system(ind1).param) ')'];
% end
for ind1=1:length(handles.system)
    temp2{ind1}=[handles.system(ind1).element ' (' num2str(handles.system(ind1).param) ')'];
end
set(handles.elementlist,'String',temp2);
% handles
guidata(hObject,handles);


% --- Executes on button press in delete_peak.
function delete_peak_Callback(hObject, eventdata, handles)
% hObject    handle to delete_peak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index_del=get(handles.elementlist,'Value');
N=length(handles.system);
if N~=0
    handles.system(index_del)=[];
    N=length(handles.system);
    if N>0
        for ind1=1:length(handles.system)
            temp{ind1}=[handles.system(ind1).element ' (' num2str(handles.system(ind1).param) ')'];
        end
        set(handles.elementlist,'String',temp);
        set(handles.elementlist,'Value',N)
    else
        set(handles.elementlist,'String','No elements added');
        set(handles.elementlist,'Value',1)
    end
else
        set(handles.elementlist,'String','No elements added');
    set(handles.elementlist,'Value',1)
end
guidata(hObject,handles);



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


% --- Executes on selection change in typeselect.
function typeselect_Callback(hObject, eventdata, handles)
% hObject    handle to typeselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns typeselect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from typeselect


% --- Executes during object creation, after setting all properties.
function typeselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to typeselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end    



% --- Executes on button press in modify.
function modify_Callback(hObject, eventdata, handles)
% hObject    handle to modify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=get(handles.elementlist,'Value');

if strcmp(handles.system(index).element,'fsp') % modify freespace
    handles.system(index).param=popup_fsp;
elseif strcmp(handles.system(index).element,'lns') % modify thin lens
    handles.system(index).param=popup_lns;
elseif strcmp(handles.system(index).element,'rff') % modify planar boundary
    handles.system(index).param=popup_rff;  
elseif strcmp(handles.system(index).element,'rfc') % modify curved boundary
    handles.system(index).param=popup_rfc;
end
N=length(handles.system);
for ind1=1:N
    temp{ind1}=[handles.system(ind1).element ' (' num2str(handles.system(ind1).param) ')'];
end
set(handles.elementlist,'String',temp);
% set(handles.elementlist,'Value',N);
guidata(hObject,handles);



% --- Executes on button press in move_up.
function move_up_Callback(hObject, eventdata, handles)
% hObject    handle to move_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=get(handles.elementlist,'Value');
if index>1
    temp2=handles.system(index-1);
    handles.system(index-1)=handles.system(index);
    handles.system(index)=temp2;
    for ind1=1:length(handles.system)
        temp{ind1}=[handles.system(ind1).element ' (' num2str(handles.system(ind1).param) ')'];
    end
    set(handles.elementlist,'String',temp);
    set(handles.elementlist,'Value',index-1);
    guidata(hObject,handles);
end

% --- Executes on button press in move_down.
function move_down_Callback(hObject, eventdata, handles)
% hObject    handle to move_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=get(handles.elementlist,'Value');
if index<length(handles.system)
    temp2=handles.system(index+1);
    handles.system(index+1)=handles.system(index);
    handles.system(index)=temp2;
    for ind1=1:length(handles.system)
        temp{ind1}=[handles.system(ind1).element ' (' num2str(handles.system(ind1).param) ')'];
    end
    set(handles.elementlist,'String',temp);
    set(handles.elementlist,'Value',index+1);
    guidata(hObject,handles);
end



function input_w0_Callback(hObject, eventdata, handles)
% hObject    handle to input_w0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_w0 as text
%        str2double(get(hObject,'String')) returns contents of input_w0 as a double
handles.input.w=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function input_w0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_w0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_R0_Callback(hObject, eventdata, handles)
% hObject    handle to input_R0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_R0 as text
%        str2double(get(hObject,'String')) returns contents of input_R0 as a double
handles.input.R=str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function input_R0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_R0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_lambda_Callback(hObject, eventdata, handles)
% hObject    handle to input_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_lambda as text
%        str2double(get(hObject,'String')) returns contents of input_lambda as a double
handles.input.lambda=str2double(get(hObject,'String'))*1e-6;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function input_lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_M2_Callback(hObject, eventdata, handles)
% hObject    handle to input_M2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_M2 as text
%        str2double(get(hObject,'String')) returns contents of input_M2 as a double
handles.input.M2=str2double(get(hObject,'String'));
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function input_M2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_M2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R_out_Callback(hObject, eventdata, handles)
% hObject    handle to R_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R_out as text
%        str2double(get(hObject,'String')) returns contents of R_out as a double
set(hObject,num2string(handles.calced.R));


% --- Executes during object creation, after setting all properties.
function R_out_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function w0_Callback(hObject, eventdata, handles)
% hObject    handle to w0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w0 as text
%        str2double(get(hObject,'String')) returns contents of w0 as a double
set(hObject,num2string(handles.calced.R));

% --- Executes during object creation, after setting all properties.
function w0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z_value_Callback(hObject, eventdata, handles)
% hObject    handle to z_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_value as text
%        str2double(get(hObject,'String')) returns contents of z_value as a double


% --- Executes during object creation, after setting all properties.
function z_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function w_out_Callback(hObject, eventdata, handles)
% hObject    handle to w_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w_out as text
%        str2double(get(hObject,'String')) returns contents of w_out as a double
set(hObject,num2string(handles.calced.w));

% --- Executes during object creation, after setting all properties.
function w_out_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lambda_out_Callback(hObject, eventdata, handles)
% hObject    handle to lambda_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda_out as text
%        str2double(get(hObject,'String')) returns contents of lambda_out as a double


% --- Executes during object creation, after setting all properties.
function lambda_out_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rayleigh_Callback(hObject, eventdata, handles)
% hObject    handle to rayleigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rayleigh as text
%        str2double(get(hObject,'String')) returns contents of rayleigh as a double


% --- Executes during object creation, after setting all properties.
function rayleigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rayleigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z_waist_Callback(hObject, eventdata, handles)
% hObject    handle to z_waist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_waist as text
%        str2double(get(hObject,'String')) returns contents of z_waist as a double


% --- Executes during object creation, after setting all properties.
function z_waist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_waist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
