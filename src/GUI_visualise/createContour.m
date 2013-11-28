function varargout = createContour(varargin)
% CREATECONTOUR MATLAB code for createContour.fig
%      CREATECONTOUR, by itself, creates a new CREATECONTOUR or raises the existing
%      singleton*.
%
%      H = CREATECONTOUR returns the handle to a new CREATECONTOUR or the handle to
%      the existing singleton*.
%
%      CREATECONTOUR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CREATECONTOUR.M with the given input arguments.
%
%      CREATECONTOUR('Property','Value',...) creates a new CREATECONTOUR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before createContour_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to createContour_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help createContour

% Last Modified by GUIDE v2.5 24-Oct-2013 14:23:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @createContour_OpeningFcn, ...
                   'gui_OutputFcn',  @createContour_OutputFcn, ...
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


% --- Executes just before createContour is made visible.
function createContour_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to createContour (see VARARGIN)

% Choose default command line output for createContour
handles.output = hObject;
handles.views = varargin{1};

%% Set up the parameters
handles.vol_axial = handles.views.axial;
handles.s_axial = size(handles.views.axial,3);
handles.current_axial = 1;
handles.h_axial = [];
handles.h_axial_position = [];
handles.vert_axial = cell(handles.s_axial,1);
handles.lines_axial = cell(handles.s_axial,1);

handles.vol_sagittal = handles.views.sagittal;
handles.s_sagittal = size(handles.views.sagittal,3); %% number of slices
handles.current_sagittal = 1;
handles.h_sagittal = [];
handles.h_sagittal_position = [];
handles.vert_sagittal = cell(handles.s_sagittal,1);
handles.lines_sagittal = cell(handles.s_sagittal,1);

handles.vol_coronal = handles.views.coronal;
handles.s_coronal = size(handles.views.coronal,3); %% number of slices
handles.current_coronal = 1;
handles.h_coronal = [];
handles.h_coronal_position = [];
handles.vert_coronal = cell(handles.s_coronal,1);
handles.lines_coronal = cell(handles.s_coronal,1);

%% Plot the data 
axes(handles.axes_logo);
imshow(imread('pictures/icon/logo.png'),[]);

axes(handles.axes_axial);
imshow(handles.views.axial(:,:,1),[]);hold on
set(handles.edit_axial,'String',['#Slices ',...
   num2str(handles.current_axial), '/',num2str(handles.s_axial)]);

axes(handles.axes_sagittal);
imshow(handles.views.sagittal(:,:,1),[]);hold on
set(handles.edit_sagittal,'String',['#Slices ',...
   num2str(handles.current_sagittal), '/',num2str(handles.s_sagittal)]);

axes(handles.axes_coronal);
imshow(handles.views.coronal(:,:,1),[]);hold on
set(handles.edit_coronal,'String',['#Slices ',...
   num2str(handles.current_coronal), '/',num2str(handles.s_coronal)]);



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes createContour wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = createContour_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all; 

% --- Executes on button press in pushbutton_ncoronal.
function pushbutton_ncoronal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ncoronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.current_coronal = handles.current_coronal + 1;

if handles.current_coronal > handles.s_coronal
    handles.current_coronal = handles.s_coronal;
end

set(handles.edit_coronal,'String',['#Slices ',...
   num2str(handles.current_coronal), '/',num2str(handles.s_coronal)]);

axes(handles.axes_coronal)
imshow(handles.views.coronal(:,:,handles.current_coronal),[]);hold on

if ~isempty(handles.vert_coronal{handles.current_coronal})

    plot(handles.vert_coronal{handles.current_coronal}(:,1),handles.vert_coronal{handles.current_coronal}(:,2),'g*','MarkerSize',10); 
    
end

axis off

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_pcoronal.
function pushbutton_pcoronal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pcoronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.current_coronal = handles.current_coronal - 1;

if handles.current_coronal < 1
    handles.current_coronal = 1;
end

set(handles.edit_coronal,'String',['#Slices ',...
   num2str(handles.current_coronal), '/',num2str(handles.s_coronal)]);

axes(handles.axes_coronal)
imshow(handles.views.coronal(:,:,handles.current_coronal),[]);hold on

if ~isempty(handles.vert_coronal{handles.current_coronal})

    plot(handles.vert_coronal{handles.current_coronal}(:,1),handles.vert_coronal{handles.current_coronal}(:,2),'g*','MarkerSize',10); 
    
end

axis off

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_nsagittal.
function pushbutton_nsagittal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_nsagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.current_sagittal = handles.current_sagittal + 1;

if handles.current_sagittal > handles.s_sagittal
    handles.current_sagittal = handles.s_sagittal;
end

set(handles.edit_sagittal,'String',['#Slices ',...
   num2str(handles.current_sagittal), '/',num2str(handles.s_sagittal)]);

axes(handles.axes_sagittal)
imshow(handles.views.sagittal(:,:,handles.current_sagittal),[]);hold on

if ~isempty(handles.vert_sagittal{handles.current_sagittal})

    plot(handles.vert_sagittal{handles.current_sagittal}(:,1),handles.vert_sagittal{handles.current_sagittal}(:,2),'g*','MarkerSize',10); 
    
end

axis off

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_psagittal.
function pushbutton_psagittal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_psagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.current_sagittal = handles.current_sagittal - 1;

if handles.current_sagittal < 1
    handles.current_sagittal = 1;
end

set(handles.edit_sagittal,'String',['#Slices ',...
   num2str(handles.current_sagittal), '/',num2str(handles.s_sagittal)]);

axes(handles.axes_sagittal)
imshow(handles.views.sagittal(:,:,handles.current_sagittal),[]);hold on

if ~isempty(handles.vert_sagittal{handles.current_sagittal})

    plot(handles.vert_sagittal{handles.current_sagittal}(:,1),handles.vert_sagittal{handles.current_sagittal}(:,2),'g*','MarkerSize',10); 
    
end

axis off

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_paxial.
function pushbutton_paxial_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_paxial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.current_axial = handles.current_axial - 1;

if handles.current_axial < 1
    handles.current_axial = 1;
end

set(handles.edit_axial,'String',['#Slices ',...
   num2str(handles.current_axial), '/',num2str(handles.s_axial)]);

axes(handles.axes_axial)
imshow(handles.views.axial(:,:,handles.current_axial),[]);hold on

if ~isempty(handles.vert_axial{handles.current_axial})

    plot(handles.vert_axial{handles.current_axial}(:,1),handles.vert_axial{handles.current_axial}(:,2),'g*','MarkerSize',10); 
    
end

axis off

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_naxial.
function pushbutton_naxial_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_naxial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.current_axial = handles.current_axial + 1;

if handles.current_axial > handles.s_axial
    handles.current_axial = handles.s_axial;
end

set(handles.edit_axial,'String',['#Slices ',...
   num2str(handles.current_axial), '/',num2str(handles.s_axial)]);

axes(handles.axes_axial)
imshow(handles.views.axial(:,:,handles.current_axial),[]);hold on

if ~isempty(handles.vert_axial{handles.current_axial})

    plot(handles.vert_axial{handles.current_axial}(:,1),handles.vert_axial{handles.current_axial}(:,2),'g*','MarkerSize',10); 
    
end


axis off

% Update handles structure
guidata(hObject, handles);


function edit_coronal_Callback(hObject, eventdata, handles)
% hObject    handle to edit_coronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_coronal as text
%        str2double(get(hObject,'String')) returns contents of edit_coronal as a double


% --- Executes during object creation, after setting all properties.
function edit_coronal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_coronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sagittal_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sagittal as text
%        str2double(get(hObject,'String')) returns contents of edit_sagittal as a double


% --- Executes during object creation, after setting all properties.
function edit_sagittal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_axial_Callback(hObject, eventdata, handles)
% hObject    handle to edit_axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_axial as text
%        str2double(get(hObject,'String')) returns contents of edit_axial as a double


% --- Executes during object creation, after setting all properties.
function edit_axial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vertices Lines] = contour_axial(handles)
verbose = 1;
nBetween = 3;

p.x = handles.h_axial_position(:,2);
p.y = handles.h_axial_position(:,1);
p.n = size(handles.h_axial_position,1);

p.x(end +1) = p.x(1);
p.y(end +1) = p.y(1);

% Interpolate to get more points
r=5;
pointsx=interp(p.x,r); pointsx=pointsx(1:end-r+1);
pointsy=interp(p.y,r); pointsy=pointsy(1:end-r+1);

axes(handles.axes_axial);
%imshow(handles.views.axial(:,:,handles.current_axial),[]);hold on
plot(pointsy,pointsx,'b','LineWidth',2); 

totalx=[]; totaly=[];
pointst = 1:length(pointsx);
i=find(pointst);
% Loop to make points evenly spaced on line pieces between landmark points
for j=1:length(i)-1
    % One line piece
    linex=pointsx(i(j):i(j+1));
    liney=pointsy(i(j):i(j+1));
    % Lenght on line through the points
    dx=[0 linex(2:end)-linex(1:end-1)];
    dy=[0 liney(2:end)-liney(1:end-1)];
    dist=cumsum(sqrt(dx.^2+dy.^2));
    % Interpolate new points evenly spaced on the line piece
    dist2=linspace(0,max(dist),nBetween);
    linex=interp1(dist,linex,dist2);
    liney=interp1(dist,liney,dist2);
    % Display the line piece
    %if(verbose),
        plot(liney,linex,'g*','MarkerSize',4);
        plot(liney(1),linex(1),'r*','MarkerSize',4);
        plot(liney(end),linex(end),'r*','MarkerSize',4); 
    %end
    % Remove Point because it is also in the next line piece
    if(j<i-1), linex(end)=[]; liney(end)=[]; end
    % Add the evenly spaced line piece to the total line
    totalx=[totalx linex];
    totaly=[totaly liney];
end
%Vertices = [totalx(:) totaly(:)];
Vertices = [p.y(:) p.x(:)];
Lines = [(1:size(Vertices,1))' ([2:size(Vertices,1) 1])'];

function [Vertices Lines] = contour_sagittal(handles)
verbose = 1;
nBetween = 3;

p.x = handles.h_sagittal_position(:,2);
p.y = handles.h_sagittal_position(:,1);
p.n = size(handles.h_sagittal_position,1);

p.x(end +1) = p.x(1);
p.y(end +1) = p.y(1);

% Interpolate to get more points
r=5;
pointsx=interp(p.x,r); pointsx=pointsx(1:end-r+1);
pointsy=interp(p.y,r); pointsy=pointsy(1:end-r+1);

axes(handles.axes_sagittal);
%imshow(handles.views.sagittal(:,:,handles.current_sagittal),[]);hold on
plot(pointsy,pointsx,'b','LineWidth',2); 

totalx=[]; totaly=[];
pointst = 1:length(pointsx);
i=find(pointst);
% Loop to make points evenly spaced on line pieces between landmark points
for j=1:length(i)-1
    % One line piece
    linex=pointsx(i(j):i(j+1));
    liney=pointsy(i(j):i(j+1));
    % Lenght on line through the points
    dx=[0 linex(2:end)-linex(1:end-1)];
    dy=[0 liney(2:end)-liney(1:end-1)];
    dist=cumsum(sqrt(dx.^2+dy.^2));
    % Interpolate new points evenly spaced on the line piece
    dist2=linspace(0,max(dist),nBetween);
    linex=interp1(dist,linex,dist2);
    liney=interp1(dist,liney,dist2);
    % Display the line piece
    %if(verbose),
        plot(liney,linex,'g*','MarkerSize',4);
        plot(liney(1),linex(1),'r*','MarkerSize',4);
        plot(liney(end),linex(end),'r*','MarkerSize',4); 
    %end
    % Remove Point because it is also in the next line piece
    if(j<i-1), linex(end)=[]; liney(end)=[]; end
    % Add the evenly spaced line piece to the total line
    totalx=[totalx linex];
    totaly=[totaly liney];
end
%Vertices = [totalx(:) totaly(:)];
Vertices = [p.y(:) p.x(:)];
Lines = [(1:size(Vertices,1))' ([2:size(Vertices,1) 1])'];

function [Vertices Lines] = contour_coronal(handles)
verbose = 1;
nBetween = 3;

p.x = handles.h_coronal_position(:,2);
p.y = handles.h_coronal_position(:,1);
p.n = size(handles.h_coronal_position,1);

p.x(end + 1) = p.x(1);
p.y(end + 1) = p.y(1);

% Interpolate to get more points
r=5;
pointsx=interp(p.x,r); pointsx=pointsx(1:end-r+1);
pointsy=interp(p.y,r); pointsy=pointsy(1:end-r+1);

axes(handles.axes_coronal);
%imshow(handles.views.coronal(:,:,handles.current_coronal),[]);hold on
plot(pointsy,pointsx,'b','LineWidth',2); 

totalx=[]; totaly=[];
pointst = 1:length(pointsx);
i=find(pointst);
% Loop to make points evenly spaced on line pieces between landmark points
for j=1:length(i)-1
    % One line piece
    linex=pointsx(i(j):i(j+1));
    liney=pointsy(i(j):i(j+1));
    % Lenght on line through the points
    dx=[0 linex(2:end)-linex(1:end-1)];
    dy=[0 liney(2:end)-liney(1:end-1)];
    dist=cumsum(sqrt(dx.^2+dy.^2));
    % Interpolate new points evenly spaced on the line piece
    dist2=linspace(0,max(dist),nBetween);
    linex=interp1(dist,linex,dist2);
    liney=interp1(dist,liney,dist2);
    % Display the line piece
    %if(verbose),
        plot(liney,linex,'g*','MarkerSize',4);
        plot(liney(1),linex(1),'r*','MarkerSize',4);
        plot(liney(end),linex(end),'r*','MarkerSize',4); 
    %end
    % Remove Point because it is also in the next line piece
    if(j<length(i)-1), linex(end)=[]; liney(end)=[]; end
    % Add the evenly spaced line piece to the total line
    totalx=[totalx linex];
    totaly=[totaly liney];
end
%Vertices = [totalx(:) totaly(:)];
Vertices = [p.y(:) p.x(:)];
Lines = [(1:size(Vertices,1))' ([2:size(Vertices,1) 1])'];



% --------------------------------------------------------------------
function uipushtool_axial_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes_axial);
handles.h_axial = impoly;
handles.h_axial_position = getPosition(handles.h_axial);
[Vertices Lines] = contour_axial(handles);

handles.vert_axial{handles.current_axial} = Vertices;
handles.lines_axial{handles.current_axial} = Lines;

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function uipushtool_sagittal_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_sagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes_sagittal);

handles.h_sagittal = impoly;
handles.h_sagittal_position = getPosition(handles.h_sagittal);
[Vertices Lines] = contour_sagittal(handles);

handles.vert_sagittal{handles.current_sagittal} = Vertices;
handles.lines_sagittal{handles.current_sagittal} = Lines;

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function uipushtool_coronal_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_coronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes_coronal);

handles.h_coronal = impoly;
handles.h_coronal_position = getPosition(handles.h_coronal);
[Vertices Lines] = contour_coronal(handles);

handles.vert_coronal{handles.current_coronal} = Vertices;
handles.lines_coronal{handles.current_coronal} = Lines;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in radiobutton_coronal.
function radiobutton_coronal_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_coronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_coronal

if (get(hObject,'Value') == get(hObject,'Max'))
    handles.vert_coronal{handles.current_coronal} = [];
    handles.lines_coronal{handles.current_coronal} = [];
	set(hObject,'Value',0);
    axes(handles.axes_coronal);
    imshow(handles.views.coronal(:,:,handles.current_coronal),[]);
else
	% Radio button is not selected-take appropriate action
end

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in radiobutton_save_coronal.
function radiobutton_save_coronal_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_save_coronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_save_coronal

if (get(hObject,'Value') == get(hObject,'Max'))
    [filename, pathname] = uiputfile('*.mat','Save as');
    if isequal(filename,0) || isequal(pathname,0)
       disp('User selected Cancel')
    else
       disp(['User selected ',fullfile(pathname,filename)])
    end
    vertex = handles.vert_coronal;
    lines = handles.lines_coronal;
    save(strcat(pathname,filename),'vertex','lines');
    set(hObject,'Value',0);
else
	% Radio button is not selected-take appropriate action
end

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in radiobutton_sagittal.
function radiobutton_sagittal_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_sagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_sagittal
if (get(hObject,'Value') == get(hObject,'Max'))
    handles.vert_sagittal{handles.current_sagittal} = [];
    handles.lines_sagittal{handles.current_sagittal} = [];
	set(hObject,'Value',0);
    axes(handles.axes_sagittal);
    imshow(handles.views.sagittal(:,:,handles.current_sagittal),[]);
else
	% Radio button is not selected-take appropriate action
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in radiobutton_save_sagittal.
function radiobutton_save_sagittal_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_save_sagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_save_sagittal
if (get(hObject,'Value') == get(hObject,'Max'))
    [filename, pathname] = uiputfile('*.mat','Save as');
    if isequal(filename,0) || isequal(pathname,0)
       disp('User selected Cancel')
    else
       disp(['User selected ',fullfile(pathname,filename)])
    end
    vertex = handles.vert_sagittal;
    lines = handles.lines_sagittal;
    save(strcat(pathname,filename),'vertex','lines');
    set(hObject,'Value',0);
else
	% Radio button is not selected-take appropriate action
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in radiobutton_axial.
function radiobutton_axial_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_axial
if (get(hObject,'Value') == get(hObject,'Max'))
    handles.vert_axial{handles.current_axial} = [];
    handles.lines_axial{handles.current_axial} = [];
	set(hObject,'Value',0);
    axes(handles.axes_axial);
    imshow(handles.views.axial(:,:,handles.current_axial),[]);
else
	% Radio button is not selected-take appropriate action
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in radiobutton_save_axial.
function radiobutton_save_axial_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_save_axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_save_axial
if (get(hObject,'Value') == get(hObject,'Max'))
    [filename, pathname] = uiputfile('*.mat','Save as');
    if isequal(filename,0) || isequal(pathname,0)
       disp('User selected Cancel')
    else
       disp(['User selected ',fullfile(pathname,filename)])
    end
    vertex = handles.vert_axial;
    lines = handles.lines_axial;
    save(strcat(pathname,filename),'vertex','lines');
    set(hObject,'Value',0);
else
	% Radio button is not selected-take appropriate action
end

% Update handles structure
guidata(hObject, handles);
