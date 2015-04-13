function varargout = main_about(varargin)
% MAIN_ABOUT MATLAB code for main_about.fig
%      MAIN_ABOUT, by itself, creates a new MAIN_ABOUT or raises the existing
%      singleton*.
%
%      H = MAIN_ABOUT returns the handle to a new MAIN_ABOUT or the handle to
%      the existing singleton*.
%
%      MAIN_ABOUT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_ABOUT.M with the given input arguments.
%
%      MAIN_ABOUT('Property','Value',...) creates a new MAIN_ABOUT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_about_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_about_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main_about

% Last Modified by GUIDE v2.5 07-Nov-2014 17:02:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_about_OpeningFcn, ...
                   'gui_OutputFcn',  @main_about_OutputFcn, ...
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


% --- Executes just before main_about is made visible.
function main_about_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main_about (see VARARGIN)

lrs_load_conf;

% Choose default command line output for main_about
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(hObject,'Name','About','NumberTitle','off');
imshow(imread(fullfile(lrs_conf.lrs_dir,'figs','lrs.png')),'parent',handles.axes_lrs);

% UIWAIT makes main_about wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main_about_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
