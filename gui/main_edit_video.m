function varargout = main_edit_video(varargin)
% MAIN_EDIT_VIDEO MATLAB code for main_edit_video.fig
%      MAIN_EDIT_VIDEO, by itself, creates a new MAIN_EDIT_VIDEO or raises the existing
%      singleton*.
%
%      H = MAIN_EDIT_VIDEO returns the handle to a new MAIN_EDIT_VIDEO or the handle to
%      the existing singleton*.
%
%      MAIN_EDIT_VIDEO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_EDIT_VIDEO.M with the given input arguments.
%
%      MAIN_EDIT_VIDEO('Property','Value',...) creates a new MAIN_EDIT_VIDEO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_edit_video_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_edit_video_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main_edit_video

% Last Modified by GUIDE v2.5 28-May-2014 22:09:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_edit_video_OpeningFcn, ...
                   'gui_OutputFcn',  @main_edit_video_OutputFcn, ...
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


% --- Executes just before main_edit_video is made visible.
function main_edit_video_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main_edit_video (see VARARGIN)

% Choose default command line output for main_edit_video
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main_edit_video wait for user response (see UIRESUME)
% uiwait(handles.figure1);

clc;
lrs_load_conf;

set(hObject,'Name','Video editor','NumberTitle','off');

set(handles.edit_input_video,'string',fullfile(lrs_conf.lrs_dir,'dataset','demo.avi'));
set(handles.edit_output_video,'string',fullfile(lrs_conf.lrs_dir,'dataset','demo_out.avi'));

img = imread(fullfile(lrs_conf.lrs_dir,'figs','no-available-image.png'));
imshow(img,'parent',handles.axes_input);


% --- Outputs from this function are returned to the command line.
function varargout = main_edit_video_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_input_video_Callback(hObject, eventdata, handles)
% hObject    handle to edit_input_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_input_video as text
%        str2double(get(hObject,'String')) returns contents of edit_input_video as a double


% --- Executes during object creation, after setting all properties.
function edit_input_video_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_input_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_show_input_video.
function pushbutton_show_input_video_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show_input_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global runningFlag;
if(runningFlag == false)
    inputVideoHandle = handles.edit_input_video;
    axesInputHandle = handles.axes_input;
    textLogMessageHandle = handles.text_log_message;
    gui_show_video(inputVideoHandle,axesInputHandle,textLogMessageHandle);
else
    disp('Stop video first!');
end


% --- Executes on button press in pushbutton_close_window.
function pushbutton_close_window_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_select_input.
function pushbutton_select_input_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lrs_load_conf;
[filename, pathname] = uigetfile({'*.avi';'*.mpg';'*.mp4';'*.*'},'File Selector',fullfile(lrs_conf.lrs_dir,'dataset'));
if(~isequal(filename,0) && ~isequal(pathname,0))
  inputVideoHandle = handles.edit_input_video;
  fullpath = fullfile(pathname, filename);
  set(inputVideoHandle,'String',fullpath);
end


% --- Executes on button press in pushbutton_stop.
function pushbutton_stop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stopFlag;
stopFlag = true;



function edit_crop_begin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_crop_begin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_crop_begin as text
%        str2double(get(hObject,'String')) returns contents of edit_crop_begin as a double


% --- Executes during object creation, after setting all properties.
function edit_crop_begin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_crop_begin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_crop_end_Callback(hObject, eventdata, handles)
% hObject    handle to edit_crop_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_crop_end as text
%        str2double(get(hObject,'String')) returns contents of edit_crop_end as a double


% --- Executes during object creation, after setting all properties.
function edit_crop_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_crop_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_resize_factor_Callback(hObject, eventdata, handles)
% hObject    handle to edit_resize_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_resize_factor as text
%        str2double(get(hObject,'String')) returns contents of edit_resize_factor as a double


% --- Executes during object creation, after setting all properties.
function edit_resize_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_resize_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_resize.
function checkbox_resize_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_resize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_resize


% --- Executes on button press in checkbox_crop.
function checkbox_crop_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_crop


% --- Executes on button press in pushbutton_process.
function pushbutton_process_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

inputVideoHandle = handles.edit_input_video;
inputFileName = get(inputVideoHandle,'String');
disp(['Input video: ' inputFileName]);

outputVideoHandle = handles.edit_output_video;
outputFileName = get(outputVideoHandle,'String');
disp(['Output video: ' outputFileName]);

if exist(inputFileName, 'file')
  displog('Loading video...');
  video_struct = load_video_file(inputFileName);
  video = convert_video_to_4d(video_struct); % show_4dvideo(video);
  video_aux = video;
  
  checkboxCropHandle = handles.checkbox_crop;
  if (get(checkboxCropHandle,'Value') == get(checkboxCropHandle,'Max'))
    displog('Cropping video...');
    
    r = str2double(get(handles.edit_crop_begin,'String'));
    s = str2double(get(handles.edit_crop_end,'String'));
    
    disp(['Begin: ', num2str(r)]);
    disp(['End: ', num2str(s)]);
    
    [video_aux] = crop_4dvideo(video,r,s);
  end

  checkboxResizeHandle = handles.checkbox_resize;
  if (get(checkboxResizeHandle,'Value') == get(checkboxResizeHandle,'Max'))
    displog('Resizing video...');
    
    f = str2double(get(handles.edit_resize_factor,'String'));
    
    disp(['Factor: ', num2str(f)]);
    
    [video_aux] = resize_4dvideo(video_aux,f);
  end

  displog('Saving video...');
  convert_video4d_to_avi(video_aux,outputFileName);
  displog('OK');
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', inputFileName);
  uiwait(msgbox(warningMessage));
end




function edit_output_video_Callback(hObject, eventdata, handles)
% hObject    handle to edit_output_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_output_video as text
%        str2double(get(hObject,'String')) returns contents of edit_output_video as a double


% --- Executes during object creation, after setting all properties.
function edit_output_video_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_output_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_show_output_video.
function pushbutton_show_output_video_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show_output_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global runningFlag;
if(runningFlag == false)
    gui_show_video(handles.edit_output_video,handles.axes_input,handles.text_log_message);
else
    disp('Stop video first!');
end
