function varargout = contours3d(varargin)
% CONTOURS3D MATLAB code for contours3d.fig
%      CONTOURS3D, by itself, creates a new CONTOURS3D or raises the existing
%      singleton*.
%
%      H = CONTOURS3D returns the handle to a new CONTOURS3D or the handle to
%      the existing singleton*.
%
%      CONTOURS3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONTOURS3D.M with the given input arguments.
%
%      CONTOURS3D('Property','Value',...) creates a new CONTOURS3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before contours3d_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to contours3d_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help contours3d

% Last Modified by GUIDE v2.5 19-Nov-2013 10:38:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @contours3d_OpeningFcn, ...
                   'gui_OutputFcn',  @contours3d_OutputFcn, ...
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


% --- Executes just before contours3d is made visible.
function contours3d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to contours3d (see VARARGIN)

%handles.orient = imread('resources/orient2.png');

% Choose default command line output for contours3d
handles.output = hObject;

handles.views = varargin{1};


handles.vol_axial = handles.views.axial;
handles.s_axial = size(handles.vol_axial,3); %% number of slices
handles.slider_axial_val = 1;%% slider value


handles.vol_sagittal = handles.views.sagittal;
handles.s_sagittal = size(handles.vol_sagittal,3); %% number of slices
handles.slider_sagittal_val = 1;%% slider value

handles.vol_coronal = handles.views.coronal;
handles.s_coronal = size(handles.vol_coronal,3); %% number of slices
handles.slider_coronal_val = 1; %% slider value

handles.view3d = 0; % flag to determine if we should show the 3dview each time the slider changes

handles.masks = cell(handles.s_axial,1); % masks, only in axial direction for the moment
handles.alpha_ch = cell(handles.s_axial,1);

axes(handles.axes_axial);
imshow(handles.vol_axial(:,:,1),[]);

set(handles.edit_AXIAL,'String',['AXIAL:     #Slices ',...
   num2str(handles.slider_axial_val), '/',num2str(handles.s_axial)]);

set(handles.slider_axial,'min',1);
set(handles.slider_axial,'max',handles.s_axial);
set(handles.slider_axial, 'SliderStep', [1/handles.s_axial , 10/handles.s_axial ]);



axes(handles.axes_sagittal);
imshow(handles.vol_sagittal(:,:,1),[]);

set(handles.edit_SAG,'String',['SAGITTAL:    #Slices ',...
   num2str(handles.slider_sagittal_val), '/',num2str(handles.s_sagittal)]);

set(handles.slider_sagittal,'min',1);
set(handles.slider_sagittal,'max',handles.s_sagittal);
set(handles.slider_sagittal, 'SliderStep', [1/handles.s_sagittal, 10/handles.s_sagittal ]);



axes(handles.axes_coronal);
imshow(handles.vol_coronal(:,:,1),[]);

set(handles.edit_COR,'String',['CORONAL:     #Slices ',...
   num2str(handles.slider_coronal_val), '/',num2str(handles.s_coronal)]);

set(handles.slider_coronal,'min',1);
set(handles.slider_coronal,'max',handles.s_coronal);
set(handles.slider_coronal, 'SliderStep', [1/handles.s_coronal , 10/handles.s_coronal ]);



axes(handles.axes_3dview);
zoom on
axis off

axes(handles.axes_orient);
imshow(imread('pictures/icon/logo.png'),[]);
 
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes contours3d wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = contours3d_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --- Executes on button press in pushbutton_3dcontours.
function pushbutton_3dcontours_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_3dcontours (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = get(handles.edit_load_mask,'String');

%% Read the images from the folder
[handles.masks handles.alpha_ch] = read_images(a,0);

% Update handles structure
guidata(hObject, handles);

% Flag 
handles.view3d = 1;

axes(handles.axes_3dview);
axis on
cla

for ax = 1%:handles.s_axial


%     %% Axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     M = zeros(4,4);
% 
%     M(1:3,4) = handles.views.axial_info{ax}.ImagePositionPatient;
%     M(4,4) = 1;
%     M(1:3,1) = handles.views.axial_info{ax}.ImageOrientationPatient(1:3).*handles.views.axial_info{ax}.PixelSpacing(1);
%     M(1:3,2) = handles.views.axial_info{ax}.ImageOrientationPatient(4:6).*handles.views.axial_info{ax}.PixelSpacing(2);
% 
%     i = [0 511];
%     j = [0 511];
% 
%     x = [];
%     y = [];
%     z = [];
% 
% 
%     for k=1:length(i)
%         for l=1:length(j)
%         p = M*[i(k) j(l) 0 1]';
% 
%             x = [x p(1)];
%             y = [y p(2)];
%             z = [z p(3)];
% 
%         end
%     end
% 
%     P1 = [x(1) y(1) z(1)];
%     P2 = [x(2) y(2) z(2)];
%     P3 = [x(3) y(3) z(3)];
%     P4 = [x(4) y(4) z(4)];
% 
%     x = [P1(1) P4(1) P2(1)];   %# [xorigin xA xB] coordinates in 3-D space
%     y = [P1(2) P4(3) P2(2)];   %# [yorigin yA yB] coordinates in 3-D space
%     z = [P1(3) P4(3) P2(3)];   %# [zorigin zA zB] coordinates in 3-D space
% 
%     origin = [512 1];  %# Vertex of triangle in image space
%     U = [0 511];       %# Vector from origin to point A in image space
%     V = [-511 511]; 
%     img = handles.masks{ax};%aa;  %# Sample image for texture map
%     
%     A = origin+U;  %# Point A
%     B = origin+V;  %# Point B
%     C = B-U;     %# Point C
% 
%     [nRows,nCols,nPages] = size(img);  %# Image dimensions
%     inputCorners = [origin; ...        %# Corner coordinates of input space
%                     A; ...
%                     B; ...
%                     C];
% 
%     outputCorners = [origin; ...      %# Corner coordinates of output space
%                      A; ...
%                      B; ...
%                      C];            
% 
%     tform = maketform('projective',...  %# Make the transformation structure
%                       inputCorners,...
%                       outputCorners);
%     triTexture = imtransform(img,tform,'bicubic',...  %# Transform the image
%                              'xdata',[1 nCols],...
%                              'ydata',[1 nRows],...
%                              'size',[nRows nCols]);
%     x = [P3(1) P4(1) P2(1) P1(1)];   %# [xorigin xA xB] coordinates in 3-D space
%     y = [P3(2) P4(2) P2(2) P1(2)];   %# [yorigin yA yB] coordinates in 3-D space
%     z = [P3(3) P4(3) P2(3) P1(3)];   %# [zorigin zA zB] coordinates in 3-D space
%     % index = [3 4; 2 1];  %# Index used to create 2-by-2 surface coordinates
%     %index = [1 4;2 3];
%     index = [4 1;3 2];
%     X1 = x(index);      %# x coordinates of surface
%     Y1 = y(index);       %# y coordinates of surface
%     Z1 = z(index);        %# z coordinates of surface


    [X,Y,Z,triTexture] = axial_rcs(handles.views.axial_info{ax},handles.masks{ax});
    hSurface = surf(X,Y,Z,'cdata',triTexture,...          %# Plot texture-mapped surface
                    'FaceColor','texturemap');hold on
                    
    set(hSurface,'edgecolor','none','facealpha','texture','alphadata',handles.alpha_ch{ax});
    alpha('direct');
    alphamap([.1;1]);

    axis equal 
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on slider movement.
function slider_coronal_Callback(hObject, eventdata, handles)
% hObject    handle to slider_coronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.slider_coronal_val = int32(get(handles.slider_coronal,'Value'));

set(handles.edit_COR,'String',['CORONAL:     #Slices ',...
   num2str(handles.slider_coronal_val), '/',num2str(handles.s_coronal)]);

axes(handles.axes_coronal);
imshow(handles.vol_coronal(:,:,handles.slider_coronal_val),[]);

% if handles.view3d == 1
%     pushbutton_3dviews_Callback(hObject, eventdata, handles);
% end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_coronal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_coronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_sagittal_Callback(hObject, eventdata, handles)
% hObject    handle to slider_sagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.slider_sagittal_val = int32(get(handles.slider_sagittal,'Value'));

set(handles.edit_SAG,'String',['SAGITTAL:    #Slices ',...
   num2str(handles.slider_sagittal_val), '/',num2str(handles.s_sagittal)]);

axes(handles.axes_sagittal);
imshow(handles.vol_sagittal(:,:,handles.slider_sagittal_val),[]);

% if handles.view3d == 1
%     pushbutton_3dviews_Callback(hObject, eventdata, handles);
% end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_sagittal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_sagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_axial_Callback(hObject, eventdata, handles)
% hObject    handle to slider_axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.slider_axial_val = int32(get(handles.slider_axial,'Value'));

set(handles.edit_AXIAL,'String',['AXIAL:     #Slices ',...
   num2str(handles.slider_axial_val), '/',num2str(handles.s_axial)]);

axes(handles.axes_axial);
imshow(handles.vol_axial(:,:,handles.slider_axial_val),[]);

% if handles.view3d == 1
%     pushbutton_3dviews_Callback(hObject, eventdata, handles);
% end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_axial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
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



function edit_AXIAL_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AXIAL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AXIAL as text
%        str2double(get(hObject,'String')) returns contents of edit_AXIAL as a double


% --- Executes during object creation, after setting all properties.
function edit_AXIAL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AXIAL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SAG_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SAG as text
%        str2double(get(hObject,'String')) returns contents of edit_SAG as a double


% --- Executes during object creation, after setting all properties.
function edit_SAG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_COR_Callback(hObject, eventdata, handles)
% hObject    handle to edit_COR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_COR as text
%        str2double(get(hObject,'String')) returns contents of edit_COR as a double


% --- Executes during object creation, after setting all properties.
function edit_COR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_COR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all; 



function edit_load_mask_Callback(hObject, eventdata, handles)
% hObject    handle to edit_load_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_load_mask as text
%        str2double(get(hObject,'String')) returns contents of edit_load_mask as a double


% --- Executes during object creation, after setting all properties.
function edit_load_mask_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_load_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
