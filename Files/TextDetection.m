clear;clc;
%% Parameter Setting
imgname = 'test7.jpg';

% Add ./image to the search path
addpath('../Images');

% Step 1 parameters setting: Parameters for Gradient Grow
MinDiversity = 0.9; % This value lies between 0 and 1
MaxVariation = 0.1; % This value lies between 0 and 1
Delta = 10;

% flag == 0 when Text is Light and Background is Dark
% flag == 1 when Text is Dark and Background is Light
flag = 0;

% Step 2 parameters setting: Parameters for Component Connection Analysis
Con_Com_small = 10; % Connected Component smallest region threshold
Con_Com_large = 2000; % Connected Component largest region threshold
Con_Com_AR = 5; % Connected Component Aspect Ratio threshold
Con_Com_EulerNumber = -3; % Connected Component Euler Number threshold

% Step 3 paramters setting: Parameters for 
stdm_thre = 0.35;

% Finally, extract the appropriate text box
T_area = 2000; % threshold in pixels


%% Step 1 Edge-Enhanced MSER
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% load the toolbox
run('./vlfeat/toolbox/vl_setup');

RGB_image = imread(imgname);
% RGB_image = imresize(RGB_image,0.4);

figure;subplot(2,3,1); imshow(RGB_image,[]);
title('Original Image');

Gray_image = uint8(rgb2gray(RGB_image));

% Extract MSER features
[r,f] = vl_mser(Gray_image,'MinDiversity',MinDiversity,'MaxVariation',MaxVariation,'Delta',Delta) ;

MSER_region = zeros(size(Gray_image));
for x=r'
    s = vl_erfill(Gray_image,x);
    MSER_region(s) = MSER_region(s)+1;
end

% figure;
% imagesc(RGB_image); 
% hold on;
% colormap gray ;
% [c,h]=contour(MSER_region,(0:max(MSER_region(:)))+.5) ;
% set(h,'color','y','linewidth',1) ;

MSER_binary = MSER_region >= 1;

% Extract Canny Edge features
Canny_binary = edge(Gray_image,'Canny');

% Combine the two features together
Canny_MSER_binary = Canny_binary & MSER_binary;

% Grow edge along gradient
[jimo,Gradient_Direction] = imgradient(Gray_image);

Growth_Direction=round((Gradient_Direction + 180) / 360 * 8); 
Growth_Direction(Growth_Direction == 0) = 8;

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
% figure; imshow(Gradient_binary);

Gra_grow_binary = ~Gradient_binary & MSER_binary;
subplot(2,3,2); imshow(Gra_grow_binary,[]);
title('Step 1: Edge-enhanced MSER');


%% Step 2 Geometric Filtering
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

Connected_comp = bwconncomp(Gra_grow_binary); % Find connected components
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

subplot(2,3,3); imshow(Masked_binary,[]);
title('Step 2: Geometric Filtering');


%% Step 3 Finding Stroke Width by Distance Transform
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

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

% Visualize the effect of stroke width filtering
subplot(2,3,4); imshow(Width_masked_binary,[]);
title('Step 3: Finding Stroke Width by Distance Transform');


%% Step 4 Morphing
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

se1=strel('disk',25);
se2=strel('disk',7);

Morphed_mask = imclose(Width_masked_binary,se1);
Morphed_mask = imopen(Morphed_mask,se2);

subplot(2,3,5); imshow(Morphed_mask,[]);
title('Step 4: Morphing');


%% Finally, extract text box
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

Connected_comp = bwconncomp(Morphed_mask);
Connected_data = regionprops(Connected_comp,'BoundingBox','Area');
boxes = round(vertcat(Connected_data(vertcat(Connected_data.Area) > T_area).BoundingBox));

subplot(2,3,6); imshow(RGB_image);
hold on;
for i=1:size(boxes,1)
    rectangle('Position', boxes(i,:),'EdgeColor','g')
end
title('Final Result!!!');
