function varargout = main(varargin)
% MAIN MATLAB code for main.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main

% Last Modified by GUIDE v2.5 07-Nov-2014 21:43:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_OpeningFcn, ...
                   'gui_OutputFcn',  @main_OutputFcn, ...
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

% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main (see VARARGIN)

% Choose default command line output for main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(hObject,'Name','LRSLibrary','NumberTitle','off');

% Close handle
set(hObject,'DeleteFcn',@gui_close);

clc;
disp('Initializing app...');

lrs_load_conf;

global method_id;
global algorithm_id;
global stopFlag;
global runningFlag;
method_id = 0;
algorithm_id = 0;
stopFlag = false;
runningFlag = false;

%img_lrs = imread('figs/lrs2.png','BackgroundColor',[0.94 0.94 0.94]);
%imshow(img_lrs,'parent',handles.axes_lrs);

imshow(imread(fullfile(lrs_conf.lrs_dir,'figs','mixed.png')),'parent',handles.axes_input);
imshow(imread(fullfile(lrs_conf.lrs_dir,'figs','outliers.png')),'parent',handles.axes_output);
imshow(imread(fullfile(lrs_conf.lrs_dir,'figs','low-rank.png')),'parent',handles.axes_lowrank);
imshow(imread(fullfile(lrs_conf.lrs_dir,'figs','sparse.png')),'parent',handles.axes_sparse);

gui_display_cputime(0,handles.axes_cputime);

set(handles.edit_input_video,'string',fullfile(lrs_conf.lrs_dir,'dataset','demo.avi'));
set(handles.edit_output_video,'string',fullfile(lrs_conf.lrs_dir,'output','output.avi'));

[list] = get_method_list();
set(handles.popupmenu_method,'string',cellstr(list));

% UIWAIT makes main wait for user response (see UIRESUME)
% uiwait(handles.figure_main);

% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_method.
function popupmenu_method_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_method
clc;
gui_display_cputime(0,handles.axes_cputime);
global method_id;
global alg_list;
contents = cellstr(get(hObject,'String'));
method = contents{get(hObject,'Value')};
set(handles.popupmenu_algorithm,'value',1);
if(method == '-')
  alg_list2{1,1} = '-';
  set(handles.popupmenu_algorithm,'string',cellstr(alg_list2));
  return;
end
disp(['Selected method: ' method]);
method_id = get_method_id(method);
disp(['method_id: ' method_id]);
alg_list = get_algorithm_list(method_id);
alg_list2 = alg_list(:,4);
set(handles.popupmenu_algorithm,'string',cellstr(alg_list2));

% --- Executes during object creation, after setting all properties.
function popupmenu_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_algorithm.
function popupmenu_algorithm_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_algorithm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_algorithm
clc;
gui_display_cputime(0,handles.axes_cputime);
global alg_list;
global algorithm_id;
contents = cellstr(get(hObject,'String'));
algorithm_desc = contents{get(hObject,'Value')}; % i.e: PCP Principal Component Pursuit (Candes et al. 2009)
if(length(algorithm_desc) > 1)
  disp(['Selected algorithm: ' algorithm_desc]);
  algorithm_info = get_algorithm_info_by_desc(algorithm_desc, alg_list); 
  algorithm_id = algorithm_info.algorithm_id; % i.e: PCP
  disp(['algorithm_id: ' algorithm_id]);
  algorithm_time = algorithm_info.algorithm_time;
  gui_display_cputime(algorithm_time,handles.axes_cputime);
end

% --- Executes on button press in pushbutton_process.
function pushbutton_process_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
global method_id;
global algorithm_id;

if(method_id == 0)
  uiwait(msgbox('Warning: select one method!'));
  return;
end

if(algorithm_id == 0)
  uiwait(msgbox('Warning: select one algorithm!'));
  return;
end

inputVideoHandle = handles.edit_input_video;
outputVideoHandle = handles.edit_output_video;

inputFileName = get(inputVideoHandle,'String');
outputFileName = get(outputVideoHandle,'String');

disp(['Input video: ' inputFileName]);
disp(['Output video: ' outputFileName]);

if exist(inputFileName, 'file')
  stats = process_video(method_id,algorithm_id,inputFileName,outputFileName);
  msg = {'Operation Completed' ...
    ' ' ['Elapsed time for decomposition: ' num2str(stats.cputime) 'sec'] ...
    ['Total elapsed time: ' num2str(stats.totaltime) 'sec'] ' '};
  uiwait(msgbox(msg,'Success'));
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', inputFileName);
  uiwait(msgbox(warningMessage));
end

% --- Executes during object creation, after setting all properties.
function popupmenu_algorithm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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


% --- Executes on button press in pushbutton_show_input.
function pushbutton_show_input_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
global runningFlag;
if(runningFlag == false)
    inputVideoHandle = handles.edit_input_video;
    axesInputHandle = handles.axes_input;
    textLogMessageHandle = handles.text_log_message;
    gui_show_video(inputVideoHandle,axesInputHandle,textLogMessageHandle);
else
    disp('Stop video first!');
end


% --- Executes on button press in pushbutton_show_output.
function pushbutton_show_output_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
global runningFlag;
if(runningFlag == false)
    outputVideoHandle = handles.edit_output_video;
    axesOutputHandle = handles.axes_output;
    textLogMessageHandle = handles.text_log_message;
    gui_show_video(outputVideoHandle,axesOutputHandle,textLogMessageHandle);
else
    disp('Stop video first!');
end

% --- Executes on button press in pushbutton_stop.
function pushbutton_stop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stopFlag;
stopFlag = true;


% --- Executes on button press in pushbutton_lowrank.
function pushbutton_lowrank_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_lowrank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
global runningFlag;
if(runningFlag == false)
  outputVideoHandle = handles.edit_output_video;
  outputFileName = get(outputVideoHandle,'String');
  L_file = gen_file_name(outputFileName,'_L');
  axesLowRankHandle = handles.axes_lowrank;
  textLogMessageHandle = handles.text_log_message;
  gui_show_video_file(L_file,axesLowRankHandle,textLogMessageHandle);
else
    disp('Stop video first!');
end

% --- Executes on button press in pushbutton_sparse.
function pushbutton_sparse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_sparse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
global runningFlag;
if(runningFlag == false)
  outputVideoHandle = handles.edit_output_video;
  outputFileName = get(outputVideoHandle,'String');
  S_file = gen_file_name(outputFileName,'_S');
  axesSparseHandle = handles.axes_sparse;
  textLogMessageHandle = handles.text_log_message;
  gui_show_video_file(S_file,axesSparseHandle,textLogMessageHandle);
else
    disp('Stop video first!');
end


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


% --- Executes on button press in pushbutton_show_results.
function pushbutton_show_results_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
lrs_load_conf;
img = imread(fullfile(lrs_conf.lrs_dir,'figs','no-available-image.png'));
imshow(img,'parent',handles.axes_input);
imshow(img,'parent',handles.axes_output);
imshow(img,'parent',handles.axes_lowrank);
imshow(img,'parent',handles.axes_sparse);
global runningFlag;
if(runningFlag == false)
  textLogMessageHandle = handles.text_log_message;  
  
  inputVideoHandle = handles.edit_input_video;
  outputVideoHandle = handles.edit_output_video;
  
  I_file = get(inputVideoHandle,'String');
  O_file = get(outputVideoHandle,'String');
  L_file = gen_file_name(O_file,'_L');
  S_file = gen_file_name(O_file,'_S');
  
  files.I_file = I_file;
  files.O_file = O_file;
  files.L_file = L_file;
  files.S_file = S_file;
  
  axes.I_axes = handles.axes_input;
  axes.O_axes = handles.axes_output;
  axes.L_axes = handles.axes_lowrank;
  axes.S_axes = handles.axes_sparse;
  
  output_extension = get_file_extension(O_file);
  if(strcmp(output_extension,'mat'))
    load(O_file);
    movs.O_mov = movobj.O;
    movs.L_mov = movobj.L;
    movs.S_mov = movobj.S;
    gui_show_results_mat(I_file, movs, axes, textLogMessageHandle);
  else
    gui_show_results(files, axes, textLogMessageHandle);
  end
else
  disp('Stop video first!');
end


% --- Executes on button press in pushbutton_gui_edit_video.
function pushbutton_gui_edit_video_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_gui_edit_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_edit_video;


% --------------------------------------------------------------------
function Main_Callback(hObject, eventdata, handles)
% hObject    handle to Main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Close_Callback(hObject, eventdata, handles)
% hObject    handle to Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;


% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function About_Callback(hObject, eventdata, handles)
% hObject    handle to About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_about;
