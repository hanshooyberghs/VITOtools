function varargout = ECMWF_Processor(varargin)
% ECMWF_Processor MATLAB code for ECMWF_Processor.fig
%      ECMWF_Processor, by itself, creates a new ECMWF_Processor or raises the existing
%      singleton*.
%
%      H = ECMWF_Processor returns the handle to a new ECMWF_Processor or the handle to
%      the existing singleton*.
%
%      ECMWF_Processor('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ECMWF_Processor.M with the given input arguments.
%
%      ECMWF_Processor('Property','Value',...) creates a new ECMWF_Processor or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ECMWF_Processor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ECMWF_Processor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ECMWF_Processor

% Last Modified by GUIDE v2.5 21-Apr-2017 16:58:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ECMWF_Processor_OpeningFcn, ...
                   'gui_OutputFcn',  @ECMWF_Processor_OutputFcn, ...
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


% --- Executes just before ECMWF_Processor is made visible.
function ECMWF_Processor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ECMWF_Processor (see VARARGIN)

% Choose default command line output for ECMWF_Processor
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

function varargout = ECMWF_Processor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    disp('reset')
    initialize_gui(gcbf, handles, true);

function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp(handles.filenames.file_2d)
disp(handles.filenames.file_3d)
if (strcmp(handles.filenames.file_2d,'') ||  strcmp(handles.filenames.file_3d,''))
    errordlg('Please set ECMWF input files','Error');
else
    set(handles.figure1, 'pointer', 'watch')
    drawnow;
    succes = 1;
    clear ME
    try 
        ecmwf2urbclim(handles.filenames,handles.coordinates,handles.data,...
        'ECMWF_data_for_urbclim.nc','log/',handles.details)
    catch ME
        succes=0;
        ME
    end
    set(handles.figure1, 'pointer', 'arrow')
    if succes ==1
        h=msgbox('Preprocessing finished','Done!');
    else
        throwException(ME)
    end
end

function surfparam_Callback(hObject, eventdata, handles)
% hObject    handle to surfparam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
[FileName,PathName] = uigetfile('*.nc','Select the NetCDF file with surface values');
handles.filenames.file_2d=[PathName,FileName];
set(handles.file2dshort,'String',FileName);
guidata(hObject,handles)

function vertprof_Callback(hObject, eventdata, handles)
% hObject    handle to vertprof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName] = uigetfile('*.nc','Select the NetCDF file with vertical profiles');
handles.filenames.file_3d=[PathName,FileName];
set(handles.file3dshort,'String',FileName);
guidata(hObject,handles)

function longitude_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
longitude = str2double(get(hObject, 'String'));
if isnan(longitude)
    set(hObject, 'String', 0);
    errordlg('Input longitude must be a number','Error');
end

if longitude> 180 || longitude < -180
    set(hObject, 'String', 0);
    errordlg('Longitude must be in range [-180,180].','Error');
    handles.coordinates.longitude = 0;
else
% Save the new volume value
handles.coordinates.longitude = longitude;
end
guidata(hObject,handles)

function longitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function latitude_Callback(hObject, eventdata, handles)
% hObject    handle to latitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of latitude as text
%        str2double(get(hObject,'String')) returns contents of latitude as a double

latitude = str2double(get(hObject, 'String'));
if isnan(latitude)
    set(hObject, 'String', 0);
    errordlg('Input latitude must be a number','Error');
end
if latitude> 90 || latitude < -90
    set(hObject, 'String', 0);
    h=errordlg('Latitude must be in range [-90,90].','Error');


else
% Save the new volume value
handles.coordinates.latitude = latitude;
end
guidata(hObject,handles)

function latitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to latitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function initialize_gui(fig_handle, handles, isreset)
 
%set coordinates
handles.coordinates.longitude = 0;
handles.coordinates.latitude  = 0;
 
set(handles.longitude, 'String', handles.coordinates.longitude);
set(handles.latitude,  'String', handles.coordinates.latitude);

% set filenames
handles.filenames.file_3d = '';
handles.filenames.file_2d = '';

set(handles.file2dshort, 'String','');
set(handles.file3dshort, 'String','');
set(handles.OverviewAdvanced, 'String','');

% provide details
handles.details.timeresolution = 1200;
handles.details.etas  = 0.472;
handles.details.psis  = 0.338;
handles.details.bch  = 6.04;
handles.details.nr_vertlevels=18;
handles.details.file_modellevels = 'modellevels.mat';

% set start and enddate
handles.data.startdate=datenum(1979,1,1);
set(handles.startdate, 'String','Start of file');

handles.data.enddate=floor(now);
set(handles.enddate, 'String','End of file');



guidata(handles.figure1, handles);

function Advanced_Callback(hObject, eventdata, handles)
% hObject    handle to Advanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

user_response = Advanced(handles.details);
hf=findobj('Name','Advanced');
close(hf)
handles.details = user_response;
set(handles.OverviewAdvanced, 'String','Advanced Data Provided');
guidata(hObject,handles)



function throwException(ME)
ME.identifier
    switch ME.identifier
        case 'InputError:InputECMWF'
            errordlg('Time format of 2d and 3d file do not match. Check input files.','Error');
        case 'InputError:CoordinatesError'
            errordlg('Coordinates do not match with input file. Check longitude and latitude.','Error');
            
        case 'InputError:TimeSelection'
            errordlg('No time steps between selected start- and enddate. Check timing.','Error');
        otherwise
            errordlg('Unknown execution error. Check input parameters.','Error');
    end


% --- Executes on button press in setstartdate.
function setstartdate_Callback(hObject, eventdata, handles)
% hObject    handle to setstartdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

user_response = guiDatePicker(handles.data.startdate);
hf=findobj('Name','guiDatePicker');
close(hf)
if user_response < handles.data.enddate
    handles.data.startdate=user_response;
    set(handles.startdate, 'String',datestr(handles.data.startdate)); 
else
    errordlg('Startdate must happen before enddate.','Error');
end
guidata(hObject,handles)


% --- Executes on button press in setenddate.
function setenddate_Callback(hObject, eventdata, handles)
% hObject    handle to setenddate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

user_response = guiDatePicker(handles.data.enddate);
hf=findobj('Name','guiDatePicker');
close(hf)
if user_response >handles.data.startdate
    handles.data.enddate=user_response;
    set(handles.enddate, 'String',datestr(handles.data.enddate));    
else
    errordlg('Startdate must happen before enddate.','Error');
end
guidata(hObject,handles)
