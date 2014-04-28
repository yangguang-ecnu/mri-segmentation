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

% for k = 1:size3
%     for i = 1:size1
%         lava_flex_sag(k,i,1:size2) = lava_flex(i,1:size2,k);
%     end
% end
% 
% for k = 1:size3
%     for j = 1:size2
%         lava_flex_cor(k,j,1:size1) = lava_flex(1:size1,j,k);
%     end
% end

for i = 1:size1
    tmp_cor(:,:,i) = imresize(lava_flex_cor(:,:,i),[size1 size2],'bicubic');
end

for i = 1:size2
    tmp_sag(:,:,i) = imresize(lava_flex_sag(:,:,i),[size1 size2],'bicubic');
end

%% Apply random transform to each view, in each plane
%% First try with rigid transforms 

% Translation
global ortho
ortho = 1;

%% First compute the transformation M, from 2D image to the 3D RCS
for i = 1:size3
    
    [lava_axM{i},lava_axM_1{i},~] = compute_M_M1(lava_flex_info{i}, ortho, 1);

end
figure
% k_ax = 1;
% aa(:,:,1) = convert2u8(lava_flex(:,:,k_ax));
% aa(:,:,2) = convert2u8(lava_flex(:,:,k_ax));
% aa(:,:,3) = convert2u8(lava_flex(:,:,k_ax));
% 
% [X, Y, Z, triTexture] = compute_RCS_v(lava_axM{k_ax}, aa);
% hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
%                'FaceColor','texturemap',...
%                'EdgeColor','none');hold on
% axis equal
% k_ax = 144;
% aa(:,:,1) = convert2u8(lava_flex(:,:,k_ax));
% aa(:,:,2) = convert2u8(lava_flex(:,:,k_ax));
% aa(:,:,3) = convert2u8(lava_flex(:,:,k_ax));
% 
% [X, Y, Z, triTexture] = compute_RCS_v(lava_axM{k_ax}, aa);
% hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
%                'FaceColor','texturemap',...
%                'EdgeColor','none');hold on
% axis equal
% k_ax = 75;
% aa(:,:,1) = convert2u8(lava_flex(:,:,k_ax));
% aa(:,:,2) = convert2u8(lava_flex(:,:,k_ax));
% aa(:,:,3) = convert2u8(lava_flex(:,:,k_ax));
% 
% [X, Y, Z, triTexture] = compute_RCS_v(lava_axM{k_ax}, aa);
% hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
%                'FaceColor','texturemap',...
%                'EdgeColor','none');hold on
% axis equal
for i = 1:size2
    
    tmp = lava_axM{1} * [size2-i 0 1]'; % size2-i+1, [i-1 0 1]
    [lava_sagM{i}, lava_sagM_1{i} ,~] = compute_M_M1_synt(lava_flex_info{1}, ortho, 2, tmp(1:3));

end
% k_sag = 215;
% ss(:,:,1) = convert2u8(lava_flex_sag(:,end:-1:1,k_sag));
% ss(:,:,2) = convert2u8(lava_flex_sag(:,end:-1:1,k_sag));
% ss(:,:,3) = convert2u8(lava_flex_sag(:,end:-1:1,k_sag));
% 
% [X, Y, Z, triTexture] = compute_RCS_v(lava_sagM{k_sag}, ss);
% hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
%                'FaceColor','texturemap',...
%                'EdgeColor','none');hold on
% axis equal
% tmp = zeros(4,1);l

for i = 1:size1
    
    tmp = lava_axM{1} * [size1-1 i-1 1]'; % [0 i-1 1]
    [lava_corM{i}, lava_corM_1{i} ,~] = compute_M_M1_synt(lava_flex_info{1}, ortho, 3, tmp(1:3));

end
% k_cor = 215;
% cc(:,:,1) = convert2u8(lava_flex_cor(:,:,k_cor));
% cc(:,:,2) = convert2u8(lava_flex_cor(:,:,k_cor));
% cc(:,:,3) = convert2u8(lava_flex_cor(:,:,k_cor));
% 
% [X, Y, Z, triTexture] = compute_RCS_v(lava_corM{k_cor}, cc);
% hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
%                'FaceColor','texturemap',...
%                'EdgeColor','none');hold on
% axis equal

%% Axial 
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
options = optimset('Display','off');

disp('--------- Calculate M^(-1) and plane eq. for each slice in the first direction -----')

figure;
for ax = 1:size(lava_flex_ax,3)
        
    [x,y,z] = calculate4corners( lava_axM{ax} );

    X_ax_v = [X_ax_v x'];
    Y_ax_v = [Y_ax_v y'];
    Z_ax_v = [Z_ax_v z'];
    
    points = [x' y' z'];
    plot3(points(:,1),points(:,2),points(:,3),'b+');hold on
    
    if ax == 1
        N1 = cross([X_ax_v(1,ax) Y_ax_v(1,ax) Z_ax_v(1,ax)]-[X_ax_v(2,ax) Y_ax_v(2,ax) Z_ax_v(2,ax)],[X_ax_v(1,ax) Y_ax_v(1,ax) Z_ax_v(1,ax)]-[X_ax_v(3,ax) Y_ax_v(3,ax) Z_ax_v(3,ax)]); % normal to the axial (ax)
        N1 = N1./norm(N1)
    end
    
end

% k_ax = 75;
% aa(:,:,1) = convert2u8(lava_flex(:,:,k_ax));
% aa(:,:,2) = convert2u8(lava_flex(:,:,k_ax));
% aa(:,:,3) = convert2u8(lava_flex(:,:,k_ax));
% 
% [X, Y, Z, triTexture] = compute_RCS_v(lava_axM{k_ax}, aa);
% hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
%                'FaceColor','texturemap',...
%                'EdgeColor','none');hold on
% axis equal

disp('--------- Calculate M^(-1) and plane eq. for each slice in the third direction -----')

for cor = 1:size1
        
    [x,y,z] = calculate4corners_cor( lava_axM{1}, lava_axM{end}, [0 cor-1;511 cor-1]); %  [cor-1 511] , [cor-1 0] [0 cor-1;511 cor-1]

    X_cor_v = [X_cor_v x'];
    Y_cor_v = [Y_cor_v y'];
    Z_cor_v = [Z_cor_v z']; 
    
    points = [x' y' z'];
    plot3(points(:,1),points(:,2),points(:,3),'r*');hold on   

    if cor == 1
        N3 = cross([X_cor_v(1,cor) Y_cor_v(1,cor) Z_cor_v(1,cor)]-[X_cor_v(2,cor) Y_cor_v(2,cor) Z_cor_v(2,cor)],[X_cor_v(1,cor) Y_cor_v(1,cor) Z_cor_v(1,cor)]-[X_cor_v(3,cor) Y_cor_v(3,cor) Z_cor_v(3,cor)]); % normal to the axial (ax)
        N3 = N3./norm(N3) 
    end
    
end

% 
% k_cor = 250;
% cc(:,:,1) = convert2u8(lava_flex_cor(:,:,k_cor));
% cc(:,:,2) = convert2u8(lava_flex_cor(:,:,k_cor));
% cc(:,:,3) = convert2u8(lava_flex_cor(:,:,k_cor));
% 
% [X, Y, Z, triTexture] = compute_RCS_img( X_cor_v(:,k_cor), Y_cor_v(:,k_cor), Z_cor_v(:,k_cor), cc);
% 
% hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
%                'FaceColor','texturemap',...
%                'EdgeColor','none');hold on
% axis equal
disp('--------- Calculate M^(-1) and plane eq. for each slice in the second direction -----')

for sag = 1:size2
        
    [x,y,z] = calculate4corners_cor( lava_axM{1}, lava_axM{end}, [511-sag 511;511-sag 0]); %  [0 511-sag;511 511-sag]

    X_sag_v = [X_sag_v x'];
    Y_sag_v = [Y_sag_v y'];
    Z_sag_v = [Z_sag_v z'];

    points = [x' y' z'];
    plot3(points(:,1),points(:,2),points(:,3),'g*');hold on   

    if sag == 1
        N2 = cross([X_sag_v(1,sag) Y_sag_v(1,sag) Z_sag_v(1,sag)]-[X_sag_v(2,sag) Y_sag_v(2,sag) Z_sag_v(2,sag)],[X_sag_v(1,sag) Y_sag_v(1,sag) Z_sag_v(1,sag)]-[X_sag_v(3,sag) Y_sag_v(3,sag) Z_sag_v(3,sag)]); % normal to the axial (ax)
        N2 = N2./norm(N2) 
    end
    
end
% 
% k_sag = 200;
% ss(:,:,1) = convert2u8(lava_flex_sag(:,:,k_sag));
% ss(:,:,2) = convert2u8(lava_flex_sag(:,:,k_sag));
% ss(:,:,3) = convert2u8(lava_flex_sag(:,:,k_sag));
% 
% [X, Y, Z, triTexture] = compute_RCS_img( X_sag_v(:,k_sag), Y_sag_v(:,k_sag), Z_sag_v(:,k_sag), ss);
% 
% hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
%                'FaceColor','texturemap',...
%                'EdgeColor','none');hold on
% axis equal

%% Select a few slices, like in T2W scans
n_slices_ax = 25;
slices_ax = 1:round(size3/n_slices_ax):size3;

st_sag   = 100;
end_sag = size2 - st_sag;
n_slices_sag = 25;
slices_sag = st_sag:ceil((end_sag-st_sag)/n_slices_sag):end_sag;

st_cor   = 150;
end_cor = size1 - st_cor;
n_slices_cor = 25;
slices_cor = st_cor:ceil((end_cor-st_cor)/n_slices_cor):end_cor;

X_ax_def = [];
Y_ax_def = [];
Z_ax_def = [];

X_sag_def = [];
Y_sag_def = [];
Z_sag_def = [];

X_cor_def = [];
Y_cor_def = [];
Z_cor_def = [];

%% Apply random deformation to the slices 
% Axial 
for i = 1:length(slices_ax)
    ind = slices_ax(i);
    out_plane = apply_t([X_ax_v(:,ind) Y_ax_v(:,ind) Z_ax_v(:,ind)],[13 12 0]');
    X_ax_def = [X_ax_def out_plane(:,1)];
    Y_ax_def = [Y_ax_def out_plane(:,2)];
    Z_ax_def = [Z_ax_def out_plane(:,3)];
%     plot3(X_ax_v(1,ind), Y_ax_v(1,ind), Z_ax_v(1,ind),'m+');hold on
%     plot3(X_ax_def(1,1), Y_ax_def(1,1), Z_ax_def(1,1), 'g*');hold on
end

% figure
% N = 9;
% k_ax = slices_ax(N);
% aa(:,:,1) = convert2u8(lava_flex(:,:,k_ax));
% aa(:,:,2) = convert2u8(lava_flex(:,:,k_ax));
% aa(:,:,3) = convert2u8(lava_flex(:,:,k_ax));
% 
% [X, Y, Z, triTexture] = compute_RCS_v(lava_axM{k_ax}, aa);
% hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
%                'FaceColor','texturemap',...
%                'EdgeColor','none');hold on
% axis equal
% 
% figure;
% aux_x = [X_ax_def(4,N);X_ax_def(3,N);X_ax_def(1,N);X_ax_def(2,N)];
% aux_y = [Y_ax_def(4,N);Y_ax_def(3,N);Y_ax_def(1,N);Y_ax_def(2,N)];
% aux_z = [Z_ax_def(4,N);Z_ax_def(3,N);Z_ax_def(1,N);Z_ax_def(2,N)];
% 
% [X, Y, Z, triTexture] = compute_RCS_img( aux_x, aux_y, aux_z, aa);
% hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
%                'FaceColor','texturemap',...
%                'EdgeColor','none');hold on
% axis equal

% Sagittal
% figure;
for i = 1:length(slices_sag)
    ind = slices_sag(i);
    out_plane = apply_t([X_sag_v(:,ind) Y_sag_v(:,ind) Z_sag_v(:,ind)],[0 13 -7]');
    X_sag_def = [X_sag_def out_plane(:,1)];
    Y_sag_def = [Y_sag_def out_plane(:,2)];
    Z_sag_def = [Z_sag_def out_plane(:,3)];
%     plot3(X_sag_v(1,ind), Y_sag_v(1,ind), Z_sag_v(1,ind),'m+');hold on
%     plot3(X_sag_def(1,1), Y_sag_def(1,1), Z_sag_def(1,1), 'g*');hold on
end
% 
% figure
% N_sag = 15;
% k_sag = slices_sag(N_sag);
% ss(:,:,1) = convert2u8(lava_flex_sag(:,:,k_sag));
% ss(:,:,2) = convert2u8(lava_flex_sag(:,:,k_sag));
% ss(:,:,3) = convert2u8(lava_flex_sag(:,:,k_sag));
% 
% [X, Y, Z, triTexture] = compute_RCS_img( X_sag_v(:,k_sag), Y_sag_v(:,k_sag), Z_sag_v(:,k_sag), ss);
% hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
%                'FaceColor','texturemap',...
%                'EdgeColor','none');hold on
% axis equal
% 
% figure
% [X, Y, Z, triTexture] = compute_RCS_img( X_sag_def(:,N), Y_sag_def(:,N), Z_sag_def(:,N), ss);
% hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
%                'FaceColor','texturemap',...
%                'EdgeColor','none');hold on
% axis equal

% Coronal
% figure;
for i = 1:length(slices_cor)
    ind = slices_cor(i);
    out_plane = apply_t([X_cor_v(:,ind) Y_cor_v(:,ind) Z_cor_v(:,ind)],[-5 0 17]');
    X_cor_def = [X_cor_def out_plane(:,1)];
    Y_cor_def = [Y_cor_def out_plane(:,2)];
    Z_cor_def = [Z_cor_def out_plane(:,3)];
%     plot3(X_sag_v(1,ind), Y_sag_v(1,ind), Z_sag_v(1,ind),'m+');hold on
%     plot3(X_sag_def(1,1), Y_sag_def(1,1), Z_sag_def(1,1), 'g*');hold on
end

% N = 9;
% k_cor = slices_cor(N);
% cc(:,:,1) = convert2u8(lava_flex_cor(:,:,k_cor));
% cc(:,:,2) = convert2u8(lava_flex_cor(:,:,k_cor));
% cc(:,:,3) = convert2u8(lava_flex_cor(:,:,k_cor));
% 
% figure;
% 
% aux_x = [X_ax_def(4,N);X_ax_def(3,N);X_ax_def(1,N);X_ax_def(2,N)];
% aux_y = [Y_ax_def(4,N);Y_ax_def(3,N);Y_ax_def(1,N);Y_ax_def(2,N)];
% aux_z = [Z_ax_def(4,N);Z_ax_def(3,N);Z_ax_def(1,N);Z_ax_def(2,N)];
% 
% [X, Y, Z, triTexture] = compute_RCS_img( aux_x, aux_y, aux_z, aa);
% hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
%                'FaceColor','texturemap',...
%                'EdgeColor','none');hold on
% axis equal
% 
% [X, Y, Z, triTexture] = compute_RCS_img( X_sag_def(:,N_sag), Y_sag_def(:,N_sag), Z_sag_def(:,N_sag), ss);
% hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
%                'FaceColor','texturemap',...
%                'EdgeColor','none');hold on
% axis equal
% 
% [X, Y, Z, triTexture] = compute_RCS_img( X_cor_def(:,N), Y_cor_def(:,N), Z_cor_def(:,N), cc);
% hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
%                'FaceColor','texturemap',...
%                'EdgeColor','none');hold on
% axis equal