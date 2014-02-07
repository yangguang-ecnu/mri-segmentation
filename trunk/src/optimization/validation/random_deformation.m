%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Read the LAVA-Flex in the 3 directions
%% and apply a random rigid transformation
%% to each one
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
serie = 7;

lava_flex      = cellimages2mat(dcmdir.dcmPatient.Study.Series(serie,1).Images);
lava_flex_info = dcmdir.dcmPatient.Study.Series(serie,1).ImagesInfo;

[size1, size2, size3] = size(lava_flex);


lava_flex_ax  = lava_flex;

for k = 1:size3
    for i = 1:size1
        lava_flex_sag(k,i,1:size2) = lava_flex(size1-i+1,size2:-1:1,k);
    end
end

for k = 1:size3
    for j = 1:size2
        lava_flex_cor(k,j,1:size1) = lava_flex(1:size1,size2-j+1,k);
    end
end

for i = 1:size1
    tmp_cor(:,:,i) = imresize(lava_flex_cor(:,:,i),[size1 size2],'bicubic');
end

for i = 1:size2
    tmp_sag(:,:,i) = imresize(lava_flex_sag(:,:,i),[size1 size2],'bicubic');
end

%% Apply random transform to each view, in each plane
%% First try with rigid transforms 

% Translation
ortho = 1;

%% First compute the transformation M, from 2D image to the 3D RCS
for i = 1:size3
    
    [lava_axM{i},lava_axM_1{i},~]     = compute_M_M1(lava_flex_info{i}, ortho, 1);
   
end

%% Axial %%
X_ax_v = [];
Y_ax_v = [];
Z_ax_v = [];

X_sag_v = [];
Y_sag_v = [];
Z_sag_v = [];

X_cor_v = [];
Y_cor_v = [];
Z_cor_v = [];


%% Axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ortho
ortho = 1;

options = optimset('Display','off');

disp('--------- Calculate M^(-1) and plane eq. for each slice in the first direction -----')
M = zeros(4,3);

for ax = 1:size(lava_flex_ax,3)
        
    [x,y,z] = calculate4corners( lava_axM{ax} );

    X_ax_v = [X_ax_v x'];
    Y_ax_v = [Y_ax_v y'];
    Z_ax_v = [Z_ax_v z'];
     
    if ax == 1
        N1 = cross([X_ax_v(1,ax) Y_ax_v(1,ax) Z_ax_v(1,ax)]-[X_ax_v(2,ax) Y_ax_v(2,ax) Z_ax_v(2,ax)],[X_ax_v(1,ax) Y_ax_v(1,ax) Z_ax_v(1,ax)]-[X_ax_v(3,ax) Y_ax_v(3,ax) Z_ax_v(3,ax)]); % normal to the axial (ax)
        N1 = N1./norm(N1)
    end
    
end

%% X_ax_v, Y_ax_v & Z_ax_v are the same for the rest of the directions,
%% but the planes are define in different way



