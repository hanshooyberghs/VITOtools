function varargout = HWDdef(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HWDdef_OpeningFcn, ...
                   'gui_OutputFcn',  @HWDdef_OutputFcn, ...
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


% --- Executes just before HWDdef is made visible.
function HWDdef_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HWDdef (see VARARGIN)

% Choose default command line output for HWDdef
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);
set(handles.figure1,'WindowStyle','modal')


% UIWAIT makes HWDdef wait for user response (see UIRESUME)
uiwait(handles.figure1);

function initialize_gui(fig_handle, handles, isreset)
% % If the metricdata field is present and the reset flag is false, it means
% % we are we are just re-initializing a GUI by calling it from the cmd line
% % while it is up. So, bail out as we dont want to reset the data.

 
handles.data.nrdays=3;
handles.data.tmin=18;
handles.data.tmax=30;
handles.data.runningmean=1;


set(handles.nrdays, 'String', handles.data.nrdays);
set(handles.tmin, 'String', handles.data.tmin);
set(handles.tmax, 'String', handles.data.tmax);
set(handles.runningmean,'value',handles.data.runningmean);

guidata(handles.figure1, handles);



% --- Outputs from this function are returned to the command line.
function varargout = HWDdef_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function tmin_Callback(hObject, eventdata, handles)
% hObject    handle to tmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tmin as text
%        str2double(get(hObject,'String')) returns contents of tmin as a double

tmin = str2double(get(hObject, 'String'));
if isnan(tmin)
    set(hObject, 'String', 0);
    errordlg('Input longitude must be a number','Error');
end

handles.data.tmin = tmin;
guidata(hObject,handles)





% --- Executes during object creation, after setting all properties.
function tmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nrdays_Callback(hObject, eventdata, handles)
% hObject    handle to nrdays (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nrdays as text
%        str2double(get(hObject,'String')) returns contents of nrdays as a double


nrdays = str2double(get(hObject, 'String'));
if isnan(nrdays)
    set(hObject, 'String', 0);
    errordlg('Input longitude must be a number','Error');
end


handles.data.nrdays = nrdays;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function nrdays_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nrdays (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tmax_Callback(hObject, eventdata, handles)
% hObject    handle to tmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tmax as text
%        str2double(get(hObject,'String')) returns contents of tmax as a double

tmax= str2double(get(hObject, 'String'));
if isnan(tmax)
    set(hObject, 'String', 0);
    errordlg('Input longitude must be a number','Error');
end


handles.data.tmax = tmax;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function tmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in confirm.
function confirm_Callback(hObject, eventdata, handles)
% hObject    handle to confirm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = handles.data;

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);


% --- Executes on button press in runningmean.
function runningmean_Callback(hObject, eventdata, handles)
% hObject    handle to runningmean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.runningmean= get(hObject,'Value');
guidata(hObject,handles)
