function varargout = UrbClimPostprocessing(varargin)
% URBCLIMPOSTPROCESSING MATLAB code for UrbClimPostprocessing.fig
%      URBCLIMPOSTPROCESSING, by itself, creates a new URBCLIMPOSTPROCESSING or raises the existing
%      singleton*.
%
%      H = URBCLIMPOSTPROCESSING returns the handle to a new URBCLIMPOSTPROCESSING or the handle to
%      the existing singleton*.
%
%      URBCLIMPOSTPROCESSING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in URBCLIMPOSTPROCESSING.M with the given input arguments.
%
%      URBCLIMPOSTPROCESSING('Property','Value',...) creates a new URBCLIMPOSTPROCESSING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UrbClimPostprocessing_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UrbClimPostprocessing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UrbClimPostprocessing

% Last Modified by GUIDE v2.5 07-May-2018 21:16:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UrbClimPostprocessing_OpeningFcn, ...
                   'gui_OutputFcn',  @UrbClimPostprocessing_OutputFcn, ...
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


% --- Executes just before UrbClimPostprocessing is made visible.
function UrbClimPostprocessing_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UrbClimPostprocessing (see VARARGIN)

% Choose default command line output for UrbClimPostprocessing
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);


function initialize_gui(fig_handle, handles, isreset)
% % If the metricdata field is present and the reset flag is false, it means
% % we are we are just re-initializing a GUI by calling it from the cmd line
% % while it is up. So, bail out as we dont want to reset the data.

% load default data
handles.indicator_data=SetUpIndicators();

% set handle data for file_in
handles.filenames.file_in = '';
set(handles.inputfile,'String','');

%set handle data for indicators
set(handles.prctile_min, 'String', handles.indicator_data.Tmin_prctile.prctile);
set(handles.prctile_mean, 'String', handles.indicator_data.Tmean_prctile.prctile);
set(handles.prctile_max, 'String', handles.indicator_data.Tmax_prctile.prctile);
set(handles.hour, 'String', handles.indicator_data.Tfixedhour_mean.hour);
set(handles.threshold_ccd, 'String', handles.indicator_data.cdd_cumul.temperature);
set(handles.threshold_tropnight, 'String', handles.indicator_data.trop_night.temperature);
set(handles.threshold_tropday, 'String', handles.indicator_data.trop_day.temperature);
set(handles.threshold_cdd_cumul_mean, 'String', handles.indicator_data.cdd_cumul_mean.temperature);
set(handles.threshold_trop_day_mean, 'String', handles.indicator_data.trop_day_mean.temperature);


%set checkbox for indicators (use default to set)
list_indicators=fieldnames(handles.indicator_data);
for i=1:length(list_indicators)
    indicator=list_indicators{i};
    checkbox=['check_',indicator];
    set(handles.(checkbox),'value',handles.indicator_data.(indicator).use);

end

% output formats
handles.output_data.tif=1;
handles.output_data.kmz=0;
handles.output_data.png=1;
%set checkbox for output formats
list_output=fieldnames(handles.output_data);
for i=1:length(list_output)
    indicator=list_output{i};
    checkbox=['check_',indicator];
    set(handles.(checkbox),'value',handles.output_data.(indicator));

end

guidata(handles.figure1, handles);


% --- Outputs from this function are returned to the command line.
function varargout = UrbClimPostprocessing_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in StartProcessing.
function StartProcessing_Callback(hObject, eventdata, handles)
% hObject    handle to StartProcessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (strcmp(handles.filenames.file_in,''))
    errordlg('Please select UrbClim output file','Error');
else
        succes=1;
        set(handles.figure1, 'pointer', 'watch')
        drawnow;
        handles.indicator_data
        try 
            ProcessUrbClim(handles.filenames.file_in,handles.indicator_data,handles.output_data)
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
        set(handles.figure1, 'pointer', 'arrow')
end

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    disp('reset')
    initialize_gui(gcbf, handles, true);


% --- Executes on button press in UrbClimOutput.
function UrbClimOutput_Callback(hObject, eventdata, handles)
% hObject    handle to UrbClimOutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName] = uigetfile('*.nc','Select the NetCDF file with UrbClim output','MultiSelect','on');
handles.filenames.file_in=[PathName,FileName];
if iscell(FileName)
    set(handles.inputfile,'String','Multiple Files');
else
    set(handles.inputfile,'String',FileName);
end
guidata(hObject,handles)



function prctile_min_Callback(hObject, eventdata, handles)
% hObject    handle to prctile_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prctile_in = str2double(get(hObject, 'String'));
if isnan(prctile_in)
    set(hObject, 'String', 0);
    errordlg('Input percentile must be a number','Error');
end
if prctile_in> 100 || prctile_in <0
    set(hObject, 'String', 0);
    h=errordlg('Percentile must be in range [0,100].','Error');

else
% Save the new volume value
handles.indicator_data.Tmin_prctile.prctile = prctile_in;
end
guidata(hObject,handles)




% --- Executes during object creation, after setting all properties.
function prctile_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prctile_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function prctile_mean_Callback(hObject, eventdata, handles)
% hObject    handle to prctile_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prctile_mean as text
%        str2double(get(hObject,'String')) returns contents of prctile_mean as a double

prctile_in = str2double(get(hObject, 'String'));
if isnan(prctile_in)
    set(hObject, 'String', 0);
    errordlg('Input percentile must be a number','Error');
end
if prctile_in> 100 || prctile_in <0
    set(hObject, 'String', 0);
    h=errordlg('Percentile must be in range [0,100].','Error');

else
% Save the new volume value
handles.indicator_data.Tmean_prctile.prctile = prctile_in;
end
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function prctile_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prctile_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function prctile_max_Callback(hObject, eventdata, handles)
% hObject    handle to prctile_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prctile_max as text
%        str2double(get(hObject,'String')) returns contents of prctile_max as a double

prctile_in = str2double(get(hObject, 'String'));
if isnan(prctile_in)
    set(hObject, 'String', 0);
    errordlg('Input percentile must be a number','Error');
end
if prctile_in> 100 || prctile_in <0
    set(hObject, 'String', 0);
    h=errordlg('Percentile must be in range [0,100].','Error');

else
% Save the new volume value
handles.indicator_data.Tmax_prctile.prctile = prctile_in;
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function prctile_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prctile_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hour_Callback(hObject, eventdata, handles)
% hObject    handle to hour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hour as text
%        str2double(get(hObject,'String')) returns contents of hour as a double

hour_in = str2double(get(hObject, 'String'));
if isnan(hour_in)
    set(hObject, 'String', 0);
    errordlg('Input percentile must be a number','Error');
end
if hour_in> 23 || hour_in <0
    set(hObject, 'String', 0);
    h=errordlg('Hour must be in range [0,23].','Error');

else
% Save the new volume value
handles.indicator_data.Tfixedhour_mean.hour = hour_in;
end
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function hour_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function threshold_ccd_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_ccd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold_ccd as text
%        str2double(get(hObject,'String')) returns contents of threshold_ccd as a double
threshold_in = str2double(get(hObject, 'String'));
if isnan(threshold_in)
    set(hObject, 'String', 0);
    errordlg('Input threshold must be a number','Error');
end
% Save the new volume value
handles.indicator_data.cdd_cumul.temperature = threshold_in;
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function threshold_ccd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold_ccd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function threshold_tropnight_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_tropnight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

threshold_in = str2double(get(hObject, 'String'));
if isnan(threshold_in)
    set(hObject, 'String', 0);
    errordlg('Input threshold must be a number','Error');
end
% Save the new volume value
handles.indicator_data.trop_night.temperature = threshold_in;
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function threshold_tropnight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold_tropnight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function threshold_tropday_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_tropday (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
threshold_in = str2double(get(hObject, 'String'));
if isnan(threshold_in)
    set(hObject, 'String', 0);
    errordlg('Input threshold must be a number','Error');
end
% Save the new volume value
handles.indicator_data.trop_day.temperature = threshold_in;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function threshold_tropday_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold_tropday (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_Tmin_mean.
function check_Tmin_mean_Callback(hObject, eventdata, handles)
% hObject    handle to check_Tmin_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.indicator_data.Tmin_mean.use= get(hObject,'Value');
guidata(hObject,handles)


% --- Executes on button press in check_Tmax_mean.
function check_Tmax_mean_Callback(hObject, eventdata, handles)
% hObject    handle to check_Tmax_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.indicator_data.Tmax_mean.use= get(hObject,'Value');
guidata(hObject,handles)


% --- Executes on button press in check_Tmean_mean.
function check_Tmean_mean_Callback(hObject, eventdata, handles)
% hObject    handle to check_Tmean_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.indicator_data.Tmean_mean.use= get(hObject,'Value');
guidata(hObject,handles)


% --- Executes on button press in check_Tmin_prctile.
function check_Tmin_prctile_Callback(hObject, eventdata, handles)
% hObject    handle to check_Tmin_prctile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.indicator_data.Tmin_prctile.use= get(hObject,'Value');
guidata(hObject,handles)


% --- Executes on button press in check_Tmean_prctile.
function check_Tmean_prctile_Callback(hObject, eventdata, handles)
% hObject    handle to check_Tmean_prctile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.indicator_data.Tmean_prctile.use= get(hObject,'Value');
guidata(hObject,handles)


% --- Executes on button press in check_Tmax_prctile.
function check_Tmax_prctile_Callback(hObject, eventdata, handles)
% hObject    handle to check_Tmax_prctile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.indicator_data.Tmax_prctile.use= get(hObject,'Value');
guidata(hObject,handles)


% --- Executes on button press in check_Tfixedhour_mean.
function check_Tfixedhour_mean_Callback(hObject, eventdata, handles)
% hObject    handle to check_Tfixedhour_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.indicator_data.Tfixedhour_mean.use= get(hObject,'Value');
guidata(hObject,handles)


% --- Executes on button press in check_cdd_cumul.
function check_cdd_cumul_Callback(hObject, eventdata, handles)
% hObject    handle to check_cdd_cumul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.indicator_data.cdd_cumul.use= get(hObject,'Value');
guidata(hObject,handles)


% --- Executes on button press in check_trop_night.
function check_trop_night_Callback(hObject, eventdata, handles)
% hObject    handle to check_trop_night (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.indicator_data.trop_night.use= get(hObject,'Value');
guidata(hObject,handles)


% --- Executes on button press in check_trop_day.
function check_trop_day_Callback(hObject, eventdata, handles)
% hObject    handle to check_trop_day (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.indicator_data.trop_day.use= get(hObject,'Value');
guidata(hObject,handles)


% --- Executes on button press in check_hwd.
function check_hwd_Callback(hObject, eventdata, handles)
% hObject    handle to check_hwd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.indicator_data.hwd.use= get(hObject,'Value');
guidata(hObject,handles)


% --- Executes on button press in SetHWDdef.
function SetHWDdef_Callback(hObject, eventdata, handles)

user_response = HWDdef();
hf=findobj('Name','HWDdef');
close(hf)
handles.indicator_data.hwd=user_response;
handles.indicator_data.hwd.use=1;
handles.indicator_data.hwd
set(handles.check_hwd,'value',handles.indicator_data.hwd.use);
guidata(hObject,handles)


% --- Executes on button press in check_png.
function check_png_Callback(hObject, eventdata, handles)
% hObject    handle to check_png (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output_data.png= get(hObject,'Value');
guidata(hObject,handles)


% --- Executes on button press in check_kmz.
function check_kmz_Callback(hObject, eventdata, handles)
% hObject    handle to check_kmz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output_data.kmz= get(hObject,'Value');
if handles.output_data.kmz == 1
    h=warndlg('Saving KMZ-files might take some time. It is not recommended to use this option for grids with more than 150 x 150 grid points.','Warning');

end
guidata(hObject,handles)


% --- Executes on button press in check_tif.
function check_tif_Callback(hObject, eventdata, handles)
% hObject    handle to check_tif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output_data.tif= get(hObject,'Value');
guidata(hObject,handles)


function throwException(ME)
ME.identifier
    switch ME.identifier
        case 'InputError:InputGeo'
            errordlg('Unmatching spatial data. Check input files.','Error');
        otherwise
            errordlg('Unknown execution error. Check input parameters.','Error');
    end



function threshold_trop_day_mean_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_tropday (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
threshold_in = str2double(get(hObject, 'String'));
if isnan(threshold_in)
    set(hObject, 'String', 0);
    errordlg('Input threshold must be a number','Error');
end
% Save the new volume value
handles.indicator_data.trop_day_mean.temperature = threshold_in;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function threshold_trop_day_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold_tropday (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function threshold_cdd_cumul_mean_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_tropday (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
threshold_in = str2double(get(hObject, 'String'));
if isnan(threshold_in)
    set(hObject, 'String', 0);
    errordlg('Input threshold must be a number','Error');
end
% Save the new volume value
handles.indicator_data.cdd_cumul_mean.temperature = threshold_in;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function threshold_cdd_cumul_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold_tropday (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in check_cdd_cumul.
function check_cdd_cumul_mean_Callback(hObject, eventdata, handles)
% hObject    handle to check_cdd_cumul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.indicator_data.cdd_cumul_mean.use= get(hObject,'Value');
guidata(hObject,handles)



% --- Executes on button press in check_trop_day.
function check_trop_day_mean_Callback(hObject, eventdata, handles)
% hObject    handle to check_trop_day (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.indicator_data.trop_day_mean.use= get(hObject,'Value');
guidata(hObject,handles)
