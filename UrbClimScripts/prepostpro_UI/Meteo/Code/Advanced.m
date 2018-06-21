function varargout = Advanced(varargin)
% ADVANCED MATLAB code for Advanced.fig
%      ADVANCED, by itself, creates a new ADVANCED or raises the existing
%      singleton*.
%
%      H = ADVANCED returns the handle to a new ADVANCED or the handle to
%      the existing singleton*.
%
%      ADVANCED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADVANCED.M with the given input arguments.
%
%      ADVANCED('Property','Value',...) creates a new ADVANCED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Advanced_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Advanced_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Advanced

% Last Modified by GUIDE v2.5 19-Apr-2017 18:00:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Advanced_OpeningFcn, ...
                   'gui_OutputFcn',  @Advanced_OutputFcn, ...
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


function Advanced_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Advanced (see VARARGIN)

% Choose default command line output for Advanced
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
handles.data = varargin{1};
initialize_gui(hObject, handles, false);
set(handles.figure1,'WindowStyle','modal')


% UIWAIT makes Advanced wait for user response (see UIRESUME)
uiwait(handles.figure1);

function varargout = Advanced_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function etas_Callback(hObject, eventdata, handles)
% hObject    handle to etas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etas as text
%        str2double(get(hObject,'String')) returns contents of etas as a double

etas = str2double(get(hObject, 'String'));
if isnan(etas)
    set(hObject, 'String', 0);
    errordlg('Input longitude must be a number','Error');
end


handles.data.etas = etas;
guidata(hObject,handles)

function etas_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function time_resolution_Callback(hObject, eventdata, handles)
% hObject    handle to time_resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_resolution as text
%        str2double(get(hObject,'String')) returns contents of time_resolution as a double


tmet = str2double(get(hObject, 'String'));
if isnan(tmet)
    set(hObject, 'String', 0);
    errordlg('Input longitude must be a number','Error');
end


handles.data.tmet = tmet;
guidata(hObject,handles)

function time_resolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function psis_Callback(hObject, eventdata, handles)
% hObject    handle to psis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of psis as text
%        str2double(get(hObject,'String')) returns contents of psis as a double

psis = str2double(get(hObject, 'String'));
if isnan(psis)
    set(hObject, 'String', 0);
    errordlg('Input longitude must be a number','Error');
end


handles.data.psis = psis;
guidata(hObject,handles)

function psis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to psis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bch_Callback(hObject, eventdata, handles)
% hObject    handle to bch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bch as text
%        str2double(get(hObject,'String')) returns contents of bch as a double

bch = str2double(get(hObject, 'String'));
if isnan(bch)
    set(hObject, 'String', 0);
    errordlg('Input longitude must be a number','Error');
end


handles.data.bch = bch;
guidata(hObject,handles)

function bch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function initialize_gui(fig_handle, handles, isreset)
 
set(handles.time_resolution, 'String', handles.data.timeresolution);
set(handles.etas, 'String', handles.data.etas);
set(handles.psis, 'String', handles.data.psis);
set(handles.bch, 'String', handles.data.bch);
set(handles.nr_vertlevels,'String',handles.data.nr_vertlevels);
set(handles.selected_modellevels,'String','Default (internal)');


handles.output=handles.data;

guidata(handles.figure1, handles);


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



function nr_vertlevels_Callback(hObject, eventdata, handles)
% hObject    handle to time_resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

input = str2double(get(hObject, 'String'));
if isnan(input)
    set(hObject, 'String', 0);
    errordlg('Input longitude must be a number','Error');
end

if input> 20
    set(hObject, 'String', 18);
    errordlg('Maximal number of levels is 20.','Error');
else
% Save the new volume value
handles.data.nr_vertlevels = input;

guidata(handles.figure1, handles);

end


% --- Executes during object creation, after setting all properties.
function nr_vertlevels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in modellevelsselect.
function modellevelsselect_Callback(hObject, eventdata, handles)
% hObject    handle to modellevelsselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName] = uigetfile('*.mat','Select the Matlab file with definition of the model levels.');
handles.data.file_modellevels=[PathName,FileName];
set(handles.selected_modellevels,'String',FileName);
guidata(hObject,handles)
