function varargout = TextDetectionDemo(varargin)
% TEXTDETECTIONDEMO MATLAB code for TextDetectionDemo.fig
%      TEXTDETECTIONDEMO, by itself, creates a new TEXTDETECTIONDEMO or raises the existing
%      singleton*.
%
%      H = TEXTDETECTIONDEMO returns the handle to a new TEXTDETECTIONDEMO or the handle to
%      the existing singleton*.
%
%      TEXTDETECTIONDEMO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEXTDETECTIONDEMO.M with the given input arguments.
%
%      TEXTDETECTIONDEMO('Property','Value',...) creates a new TEXTDETECTIONDEMO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TextDetectionDemo_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TextDetectionDemo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TextDetectionDemo

% Last Modified by GUIDE v2.5 16-May-2014 10:22:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TextDetectionDemo_OpeningFcn, ...
                   'gui_OutputFcn',  @TextDetectionDemo_OutputFcn, ...
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


% --- Executes just before TextDetectionDemo is made visible.
function TextDetectionDemo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TextDetectionDemo (see VARARGIN)

% Choose default command line output for TextDetectionDemo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TextDetectionDemo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TextDetectionDemo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename ,pathname]=uigetfile({'*.jpg';'*.bmp';'*.tif';'*.gif'},'Choose Image');
if filename==0
    return;
else
    str=[pathname filename];
    global RGB_image;
    RGB_image=imread(str);
    axes(handles.axes1);
    cla reset;
    imshow(RGB_image,[]);
end

% --------------------------------------------------------------------
function HelpMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to HelpMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox({'This Demo is developed by Yinjie Huang from University of Central Florida Machine Learning Lab.';'';...
    'This Demo is used to show Automatic Text Detection.';'';...
    'For more information, please refer to the ReadMe file!!'});

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RGB_image;
run('./Files/vlfeat/toolbox/vl_setup');
Gray_image = uint8(rgb2gray(RGB_image));
MinDiversity = 0.9; % This value lies between 0 and 1
MaxVariation = 0.1; % This value lies between 0 and 1
Delta = 10;
% Extract MSER features
[r,f] = vl_mser(Gray_image,'MinDiversity',MinDiversity,'MaxVariation',MaxVariation,'Delta',Delta) ;

MSER_region = zeros(size(Gray_image));
for x=r'
    s = vl_erfill(Gray_image,x);
    MSER_region(s) = MSER_region(s)+1;
end
MSER_binary = MSER_region >= 1;

% Extract Canny Edge features
Canny_binary = edge(Gray_image,'Canny');

% Combine the two features together
Canny_MSER_binary = Canny_binary & MSER_binary;

% Grow edge along gradient
[jimo,Gradient_Direction] = imgradient(Gray_image);

Growth_Direction=round((Gradient_Direction + 180) / 360 * 8); 
Growth_Direction(Growth_Direction == 0) = 8;

flag=str2num(get(handles.edit1,'String'));
if flag==1 % when we have dark text and light background, reverse growth direction
    Growth_Direction=mod(Growth_Direction + 3, 8) + 1;
end

% Bulid growing structure template
% Diagonal template
NW_Template = [1,1,1,1,1,0,0;1,1,1,1,1,0,0;1,1,1,1,1,0,0;...
    1,1,1,1,0,0,0;1,1,1,0,0,0,0;zeros(2,7)];

% Horizontal and Vertical template
N_Template = [0,1,1,1,1,1,0;0,1,1,1,1,1,0;0,1,1,1,1,1,0;...
    0,1,1,1,1,1,0;zeros(3,7)];

N = strel(N_Template);
W = strel(rot90(N_Template,1));
S = strel(rot90(N_Template,2));
E = strel(rot90(N_Template,3));
NW = strel(NW_Template);
SW = strel(rot90(NW_Template,1));
SE = strel(rot90(NW_Template,2));
NE = strel(rot90(NW_Template,3));

Strels = [NE,N,NW,W,SW,S,SE,E];

Gradient_binary = false(size(Canny_MSER_binary));
% Let's grow edges
for i = 1:numel(Strels)
    Temp_Direction = false(size(Canny_MSER_binary));
    Temp_Direction(Canny_MSER_binary == true & Growth_Direction == i ) = true;
    Temp_Direction = imdilate(Temp_Direction,Strels(i));
    Gradient_binary = Gradient_binary | Temp_Direction;
end
global Gra_grow_binary;
Gra_grow_binary = ~Gradient_binary & MSER_binary;
axes(handles.axes2);
cla reset;
imshow(Gra_grow_binary,[]);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Gra_grow_binary;
Connected_comp = bwconncomp(Gra_grow_binary); % Find connected components
Con_Com_small = 10; % Connected Component smallest region threshold
Con_Com_large = 2000; % Connected Component largest region threshold
Con_Com_AR = 5; % Connected Component Aspect Ratio threshold
Con_Com_EulerNumber = -3; % Connected Component Euler Number threshold
Connected_data = regionprops(Connected_comp,'Area','EulerNumber','Image');
global Masked_binary;
Masked_binary = Gra_grow_binary;

% Rule 1, reject very large or very small objects
Masked_binary(vertcat(Connected_comp.PixelIdxList{[Connected_data.Area] < Con_Com_small | [Connected_data.Area] > Con_Com_large})) = 0;
% Rule 2, reject very large or very small Aspect Ratio
C = {Connected_data(:).Image};
AspectRatio = cellfun(@(x)max(size(x))./min(size(x)),C);
Masked_binary(vertcat(Connected_comp.PixelIdxList{AspectRatio >= Con_Com_AR})) = 0;
% Rule 3, reject objects with lots of holes
Masked_binary(vertcat(Connected_comp.PixelIdxList{[Connected_data.EulerNumber] <= Con_Com_EulerNumber})) = 0;

axes(handles.axes3);
cla reset;
imshow(Masked_binary,[]);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Masked_binary;
stdm_thre=str2double(get(handles.edit2,'String'));
Distance_image = bwdist(~Masked_binary);

% The following algorithm follows the paper
Distance_image = round(Distance_image);

% Define 8-connected neighbors
connect = [1 0;-1 0;1 1;0 1;-1 1;1 -1;0 -1;-1 -1]';
padded_distance_image = padarray(Distance_image,[1,1]);
D_ind = find(padded_distance_image ~= 0);
sz=size(padded_distance_image);

% Compare current pixel with its neighbors
neighbor_ind = repmat(D_ind,[1,8]);
[x,y] = ind2sub(sz,neighbor_ind);
x = bsxfun(@plus,x,connect(1,:));
y = bsxfun(@plus,y,connect(2,:));
neighbor_ind = sub2ind(sz,x,y);
LookUp = bsxfun(@lt,padded_distance_image(neighbor_ind),padded_distance_image(D_ind));

% Propagate local maximum stroke values to neighbors recursively
MaxStroke = max(max(padded_distance_image));
for Stroke = MaxStroke:-1:1
    neighbor_ind_temp = neighbor_ind(padded_distance_image(D_ind) == Stroke,:);
    LookUp_temp = LookUp(padded_distance_image(D_ind) == Stroke,:);
    NeighborIndex = neighbor_ind_temp(LookUp_temp);
    while ~isempty(NeighborIndex)
        padded_distance_image(NeighborIndex) = Stroke;       
        [~,ia,~] = intersect(D_ind,NeighborIndex);
        neighbor_ind_temp = neighbor_ind(ia,:);
        LookUp_temp = LookUp(ia,:);
        NeighborIndex = neighbor_ind_temp(LookUp_temp);
    end
end

% Remove pad to restore original image size
Width_image = padded_distance_image(2:end-1,2:end-1);

% figure; imshow(Width_image);
% caxis([0 max(max(Width_image))]); axis image, colormap('jet'), colorbar;

Connected_comp = bwconncomp(Masked_binary);
global Width_masked_binary;
Width_masked_binary = Masked_binary;
for i = 1:Connected_comp.NumObjects
    strokewidths = Width_image(Connected_comp.PixelIdxList{i});
    % Compute normalized stroke width variation and compare to common value
    if std(strokewidths)/mean(strokewidths) > stdm_thre
        Width_masked_binary(Connected_comp.PixelIdxList{i}) = 0; % Remove from text candidates
    end
end
axes(handles.axes4);
cla reset;
imshow(Width_masked_binary,[]);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Width_masked_binary;
se1=strel('disk',25);
se2=strel('disk',7);
global Morphed_mask;

Morphed_mask = imclose(Width_masked_binary,se1);
Morphed_mask = imopen(Morphed_mask,se2);
axes(handles.axes5);
cla reset;
imshow(Morphed_mask,[]);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Morphed_mask;
T_area=str2num(get(handles.edit3,'String'));

Connected_comp = bwconncomp(Morphed_mask);
Connected_data = regionprops(Connected_comp,'BoundingBox','Area');
boxes = round(vertcat(Connected_data(vertcat(Connected_data.Area) > T_area).BoundingBox));

axes(handles.axes12);
cla reset;
global RGB_image;
imshow(RGB_image,[]);
hold on;
for i=1:size(boxes,1)
    rectangle('Position', boxes(i,:),'EdgeColor','g')
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RGB_image;
run('./Files/vlfeat/toolbox/vl_setup');
Gray_image = uint8(rgb2gray(RGB_image));
MinDiversity = 0.9; % This value lies between 0 and 1
MaxVariation = 0.1; % This value lies between 0 and 1
Delta = 10;
% Extract MSER features
[r,f] = vl_mser(Gray_image,'MinDiversity',MinDiversity,'MaxVariation',MaxVariation,'Delta',Delta) ;

MSER_region = zeros(size(Gray_image));
for x=r'
    s = vl_erfill(Gray_image,x);
    MSER_region(s) = MSER_region(s)+1;
end
MSER_binary = MSER_region >= 1;

% Extract Canny Edge features
Canny_binary = edge(Gray_image,'Canny');

% Combine the two features together
Canny_MSER_binary = Canny_binary & MSER_binary;

% Grow edge along gradient
[jimo,Gradient_Direction] = imgradient(Gray_image);

Growth_Direction=round((Gradient_Direction + 180) / 360 * 8); 
Growth_Direction(Growth_Direction == 0) = 8;

flag=str2num(get(handles.edit1,'String'));
if flag==1 % when we have dark text and light background, reverse growth direction
    Growth_Direction=mod(Growth_Direction + 3, 8) + 1;
end

% Bulid growing structure template
% Diagonal template
NW_Template = [1,1,1,1,1,0,0;1,1,1,1,1,0,0;1,1,1,1,1,0,0;...
    1,1,1,1,0,0,0;1,1,1,0,0,0,0;zeros(2,7)];

% Horizontal and Vertical template
N_Template = [0,1,1,1,1,1,0;0,1,1,1,1,1,0;0,1,1,1,1,1,0;...
    0,1,1,1,1,1,0;zeros(3,7)];

N = strel(N_Template);
W = strel(rot90(N_Template,1));
S = strel(rot90(N_Template,2));
E = strel(rot90(N_Template,3));
NW = strel(NW_Template);
SW = strel(rot90(NW_Template,1));
SE = strel(rot90(NW_Template,2));
NE = strel(rot90(NW_Template,3));

Strels = [NE,N,NW,W,SW,S,SE,E];

Gradient_binary = false(size(Canny_MSER_binary));
% Let's grow edges
for i = 1:numel(Strels)
    Temp_Direction = false(size(Canny_MSER_binary));
    Temp_Direction(Canny_MSER_binary == true & Growth_Direction == i ) = true;
    Temp_Direction = imdilate(Temp_Direction,Strels(i));
    Gradient_binary = Gradient_binary | Temp_Direction;
end
Gra_grow_binary = ~Gradient_binary & MSER_binary;
axes(handles.axes2);
cla reset;
imshow(Gra_grow_binary,[]);

Connected_comp = bwconncomp(Gra_grow_binary); % Find connected components
Con_Com_small = 10; % Connected Component smallest region threshold
Con_Com_large = 2000; % Connected Component largest region threshold
Con_Com_AR = 5; % Connected Component Aspect Ratio threshold
Con_Com_EulerNumber = -3; % Connected Component Euler Number threshold
Connected_data = regionprops(Connected_comp,'Area','EulerNumber','Image');
Masked_binary = Gra_grow_binary;

% Rule 1, reject very large or very small objects
Masked_binary(vertcat(Connected_comp.PixelIdxList{[Connected_data.Area] < Con_Com_small | [Connected_data.Area] > Con_Com_large})) = 0;
% Rule 2, reject very large or very small Aspect Ratio
C = {Connected_data(:).Image};
AspectRatio = cellfun(@(x)max(size(x))./min(size(x)),C);
Masked_binary(vertcat(Connected_comp.PixelIdxList{AspectRatio >= Con_Com_AR})) = 0;
% Rule 3, reject objects with lots of holes
Masked_binary(vertcat(Connected_comp.PixelIdxList{[Connected_data.EulerNumber] <= Con_Com_EulerNumber})) = 0;

axes(handles.axes3);
cla reset;
imshow(Masked_binary,[]);

stdm_thre=str2double(get(handles.edit2,'String'));
Distance_image = bwdist(~Masked_binary);

% The following algorithm follows the paper
Distance_image = round(Distance_image);

% Define 8-connected neighbors
connect = [1 0;-1 0;1 1;0 1;-1 1;1 -1;0 -1;-1 -1]';
padded_distance_image = padarray(Distance_image,[1,1]);
D_ind = find(padded_distance_image ~= 0);
sz=size(padded_distance_image);

% Compare current pixel with its neighbors
neighbor_ind = repmat(D_ind,[1,8]);
[x,y] = ind2sub(sz,neighbor_ind);
x = bsxfun(@plus,x,connect(1,:));
y = bsxfun(@plus,y,connect(2,:));
neighbor_ind = sub2ind(sz,x,y);
LookUp = bsxfun(@lt,padded_distance_image(neighbor_ind),padded_distance_image(D_ind));

% Propagate local maximum stroke values to neighbors recursively
MaxStroke = max(max(padded_distance_image));
for Stroke = MaxStroke:-1:1
    neighbor_ind_temp = neighbor_ind(padded_distance_image(D_ind) == Stroke,:);
    LookUp_temp = LookUp(padded_distance_image(D_ind) == Stroke,:);
    NeighborIndex = neighbor_ind_temp(LookUp_temp);
    while ~isempty(NeighborIndex)
        padded_distance_image(NeighborIndex) = Stroke;       
        [~,ia,~] = intersect(D_ind,NeighborIndex);
        neighbor_ind_temp = neighbor_ind(ia,:);
        LookUp_temp = LookUp(ia,:);
        NeighborIndex = neighbor_ind_temp(LookUp_temp);
    end
end

% Remove pad to restore original image size
Width_image = padded_distance_image(2:end-1,2:end-1);

% figure; imshow(Width_image);
% caxis([0 max(max(Width_image))]); axis image, colormap('jet'), colorbar;

Connected_comp = bwconncomp(Masked_binary);
Width_masked_binary = Masked_binary;
for i = 1:Connected_comp.NumObjects
    strokewidths = Width_image(Connected_comp.PixelIdxList{i});
    % Compute normalized stroke width variation and compare to common value
    if std(strokewidths)/mean(strokewidths) > stdm_thre
        Width_masked_binary(Connected_comp.PixelIdxList{i}) = 0; % Remove from text candidates
    end
end
axes(handles.axes4);
cla reset;
imshow(Width_masked_binary,[]);

se1=strel('disk',25);
se2=strel('disk',7);

Morphed_mask = imclose(Width_masked_binary,se1);
Morphed_mask = imopen(Morphed_mask,se2);
axes(handles.axes5);
cla reset;
imshow(Morphed_mask,[]);

T_area=str2num(get(handles.edit3,'String'));

Connected_comp = bwconncomp(Morphed_mask);
Connected_data = regionprops(Connected_comp,'BoundingBox','Area');
boxes = round(vertcat(Connected_data(vertcat(Connected_data.Area) > T_area).BoundingBox));

axes(handles.axes12);
cla reset;
global RGB_image;
imshow(RGB_image,[]);
hold on;
for i=1:size(boxes,1)
    rectangle('Position', boxes(i,:),'EdgeColor','g')
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
