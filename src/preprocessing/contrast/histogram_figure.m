function varargout = histogram_figure(varargin)
% HISTOGRAM_FIGURE M-file for histogram_figure.fig
%      HISTOGRAM_FIGURE, by itself, creates a new HISTOGRAM_FIGURE or raises the existing
%      singleton*.
%
%      H = HISTOGRAM_FIGURE returns the handle to a new HISTOGRAM_FIGURE or the handle to
%      the existing singleton*.
%
%      HISTOGRAM_FIGURE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HISTOGRAM_FIGURE.M with the given input arguments.
%
%      HISTOGRAM_FIGURE('Property','Value',...) creates a new HISTOGRAM_FIGURE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before histogram_figure_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to histogram_figure_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help histogram_figure

% Last Modified by GUIDE v2.5 26-Apr-2011 11:44:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @histogram_figure_OpeningFcn, ...
                   'gui_OutputFcn',  @histogram_figure_OutputFcn, ...
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


% --- Executes just before histogram_figure is made visible.
function histogram_figure_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to histogram_figure (see VARARGIN)

% Choose default command line output for histogram_figure
handles.output = hObject;

axes(handles.IMAGE_AXES);
axis off

axes(handles.IMAGE_EQ_AXES);
axis off

% axes(handles.HISTOGRAM_AXES);
% axis off
% 
% axes(handles.HISTOGRAM_EQ_AXES);
% axis off
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes histogram_figure wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = histogram_figure_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in BARRAS_PUSHBUTTON.
function BARRAS_PUSHBUTTON_Callback(hObject, eventdata, handles)
% hObject    handle to BARRAS_PUSHBUTTON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Plot the histogram
x = getappdata(histogram_figure, 'Histogram_x');
h = getappdata(histogram_figure, 'Histogram');
n = getappdata(histogram_figure, 'Image_channels');
plot_axes_bar(handles.HISTOGRAM_AXES,x,h,n);

% Plot the equalize histogram
h_eq = getappdata(histogram_figure, 'Histogram_eq');
plot_axes_bar(handles.HISTOGRAM_EQ_AXES,x,h_eq,n);



% --- Executes on button press in STEM_PUSHBUTTON.
function STEM_PUSHBUTTON_Callback(hObject, eventdata, handles)
% hObject    handle to STEM_PUSHBUTTON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Plot the histogram
x = getappdata(histogram_figure, 'Histogram_x');
h = getappdata(histogram_figure, 'Histogram');
n = getappdata(histogram_figure, 'Image_channels');
plot_axes_st(handles.HISTOGRAM_AXES,x,h,n);

% Plot the equalize histogram
h_eq = getappdata(histogram_figure, 'Histogram_eq');
plot_axes_st(handles.HISTOGRAM_EQ_AXES,x,h_eq,n);

% --- Executes on button press in PLOT_PUSHBUTTON.
function PLOT_PUSHBUTTON_Callback(hObject, eventdata, handles)
% hObject    handle to PLOT_PUSHBUTTON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Plot the histogram
x = getappdata(histogram_figure, 'Histogram_x');
h = getappdata(histogram_figure, 'Histogram');
n = getappdata(histogram_figure, 'Image_channels');
plot_axes_p(handles.HISTOGRAM_AXES,x,h,n);

% Plot the equalize histogram
h_eq = getappdata(histogram_figure, 'Histogram_eq');
plot_axes_p(handles.HISTOGRAM_EQ_AXES,x,h_eq,n);


% --- Executes on button press in HISTOGRAM_PUSHBUTTOM.
function HISTOGRAM_PUSHBUTTOM_Callback(hObject, eventdata, handles)
% hObject    handle to HISTOGRAM_PUSHBUTTOM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = getappdata(histogram_figure, 'Image');
[rows cols] = size(I);
[h,x,n] = histogram_image(I);
setappdata(histogram_figure, 'Image_channels', n);
setappdata(histogram_figure, 'Histogram', h);
setappdata(histogram_figure, 'Histogram_x', x);
plot_axes(handles.HISTOGRAM_AXES,x,h./(rows*cols),n);
setappdata(histogram_figure, 'Axes_Histogram', handles.HISTOGRAM_AXES);

% --- Executes on button press in HISTOGRAM_EQ_PUSHBUTTON.
function HISTOGRAM_EQ_PUSHBUTTON_Callback(hObject, eventdata, handles)
% hObject    handle to HISTOGRAM_EQ_PUSHBUTTON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = getappdata(histogram_figure, 'Image');
[rows cols] = size(I);
h = getappdata(histogram_figure, 'Histogram');
[I_eq,h_eq,cdf,h_]=histogram_eq_image(I,h);
% Save the inputs as equalize image & equalize histogram 
setappdata(histogram_figure, 'Image_eq', I_eq);
setappdata(histogram_figure, 'Histogram_eq', h_eq);
setappdata(histogram_figure, 'Cdf', cdf);
setappdata(histogram_figure, 'Histogram_1', h_);
% Plot the equalize histogram
x = getappdata(histogram_figure, 'Histogram_x');
n = getappdata(histogram_figure, 'Image_channels');

plot_axes(handles.HISTOGRAM_EQ_AXES,x,h_eq./(rows*cols),n);
setappdata(histogram_figure, 'Axes_Histogram_EQ', handles.HISTOGRAM_EQ_AXES);

% --- Executes on button press in IMAGE_EQ_PUSHBUTTON.
function IMAGE_EQ_PUSHBUTTON_Callback(hObject, eventdata, handles)
% hObject    handle to IMAGE_EQ_PUSHBUTTON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I_eq = getappdata(histogram_figure, 'Image_eq');
axes(handles.IMAGE_EQ_AXES)
imshow(I_eq,[])

% --- Executes on button press in SELECTIMAGE_PUSHBUTTON.
function SELECTIMAGE_PUSHBUTTON_Callback(hObject, eventdata, handles)
% hObject    handle to SELECTIMAGE_PUSHBUTTON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName, FilePath] = uigetfile({ '*.jpg'; '*.jpeg'; '*.bmp'; '*.gif'; '*.tiff'; '*.png' }, 'Select an Image...');
% read the image
I = imread([FilePath, FileName]);
% set the image in the user data, if this is not done, the image would be
% lost
setappdata(histogram_figure, 'Image', I);

axes(handles.IMAGE_AXES)
imshow(I,[])
axis off


% --- Executes on button press in PDF_CHECKBOX.
function PDF_CHECKBOX_Callback(hObject, eventdata, handles)
% hObject    handle to PDF_CHECKBOX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PDF_CHECKBOX


% --- Executes on button press in CDF_CHECKBOX.
function CDF_CHECKBOX_Callback(hObject, eventdata, handles)
% hObject    handle to CDF_CHECKBOX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cdf = getappdata(histogram_figure, 'Cdf');
h1 = getappdata(histogram_figure, 'Histogram_1');
x = getappdata(histogram_figure, 'Histogram_x');
I = getappdata(histogram_figure, 'Image');

[rows cols] = size(I);
% Hint: get(hObject,'Value') returns toggle state of CDF_CHECKBOX
if (get(hObject,'Value') == get(hObject,'Max'))
   % Checkbox is checked-take appropriate action
    axes(handles.HISTOGRAM_AXES);
    hold on
    plot_cdf = plot(x,cdf./(rows*cols),'k');
    set(plot_cdf,'Visible','on');
    %plot_cdf = plot(x,cdf./(rows*cols),'k');
    %setappdata(histogram_figure, 'Plot_cdf', plot_cdf);
    %axis([0 max(x) 0 1])
    axes(handles.HISTOGRAM_EQ_AXES);
    hold on
    plot_h1 = plot(x,h1./(rows*cols),'k');
    set(plot_h1,'Visible','on');
    %plot_h1 = plot(x,h1./(rows*cols),'k');
    %setappdata(histogram_figure, 'Plot_h1', plot_h1);
    %axis([0 max(x) 0 1])
else
    %plot_cdf = getappdata(histogram_figure, 'Plot_cdf');
    %plot_h1 = getappdata(histogram_figure, 'Plot_h1');
    %axes(handles.HISTOGRAM_AXES);
    %hold on
    %plot_cdf = plot(x,cdf./(rows*cols),'k');
    axes(handles.HISTOGRAM_EQ_AXES);
    plot_h1 = plot(x,h1./(rows*cols),'k');
    %set(plot_cdf,'Visible','off');
    %delete plot_cdf
    %axes(handles.HISTOGRAM_EQ_AXES);
    %hold on
    %plot_h1 = plot(x,h1./(rows*cols),'k');
    %set(plot_h1,'Visible','off');
    axes(handles.HISTOGRAM_EQ_AXES);
    plot_h1 = plot(x,h1./(rows*cols),'k');
   % Checkbox is not checked-take appropriate action
end

% --- Executes on selection change in RGB_POPUPMENU.
function RGB_POPUPMENU_Callback(hObject, eventdata, handles)
% hObject    handle to RGB_POPUPMENU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns RGB_POPUPMENU contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RGB_POPUPMENU
% Determine the selected data set.
str = get(hObject, 'String');
val = get(hObject,'Value');
n = getappdata(histogram_figure, 'Image_channels');

if n == 3
    x = getappdata(histogram_figure, 'Histogram_x');
    h = getappdata(histogram_figure, 'Histogram');
    h_eq = getappdata(histogram_figure, 'Histogram_eq');
    % Set current data to the selected data set.
    switch str{val};
        case 'RGB-channel' % User selects peaks.
            % Plot the histogram
            plot_axes(handles.HISTOGRAM_AXES,x,h,n);
            % Plot the equalize histogram
            plot_axes(handles.HISTOGRAM_EQ_AXES,x,h_eq,n);

        case 'red'
            % Plot the histogram
            plot_axes(handles.HISTOGRAM_AXES,x,h(:,1),1,'r');
            % Plot the equalize histogram
            plot_axes(handles.HISTOGRAM_EQ_AXES,x,h_eq(:,1),1,'r');

        case 'green'
            % Plot the histogram
            plot_axes(handles.HISTOGRAM_AXES,x,h(:,2),1,'g');
            % Plot the equalize histogram
            plot_axes(handles.HISTOGRAM_EQ_AXES,x,h_eq(:,2),1,'g');

        case 'blue'
            % Plot the histogram
            plot_axes(handles.HISTOGRAM_AXES,x,h(:,3),1,'b');
            % Plot the equalize histogram
            plot_axes(handles.HISTOGRAM_EQ_AXES,x,h_eq(:,3),1,'b');

    end
end
            

% --- Executes during object creation, after setting all properties.
function RGB_POPUPMENU_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RGB_POPUPMENU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
