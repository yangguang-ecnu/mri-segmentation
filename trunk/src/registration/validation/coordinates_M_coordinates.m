function st = coordinates_M_coordinates(lava_flex_n, lava_flex_info)

%% Generate the deformed MRI

[size1, size2, size3] = size(lava_flex_n);

disp('--------- Define each direction -----')
lava_flex_ax  = lava_flex_n;

for k = 1:size3
    for i = 1:size1
        lava_flex_sag(k,i,:) = lava_flex_n(i,:,k); %size1-i+1,size2:-1:1,k
        %         lava_flex_sag(k,i,1:size2) = lava_flex(i,1:size2,k); %size1-i+1,size2:-1:1,k
    end
end


for k = 1:size3
    for j = 1:size2
        lava_flex_cor(k,j,1:size1) = lava_flex_n(1:size1,j,k); % 1:size1 size2-j+1
    end
end

%% Apply random transform to each view, in each plane
%% First try with rigid transforms

% Translation

ortho = 1;

disp('--------- Calculate transformations -----')
%% First compute the transformation M, from 2D image to the 3D RCS
for i = 1:size3
    
    [lava_axM{i},lava_axM_1{i},~] = compute_M_M1(lava_flex_info{i}, ortho, 1);
    
end


for i = 1:size2
    
    tmp = lava_axM{1} * [i-1 0 1]';
    [lava_sagM{i}, lava_sagM_1{i} ,~] = compute_M_M1_synt(lava_flex_info{1}, ortho, 2, tmp(1:3));
    
end

for i = 1:size1
    
    tmp = lava_axM{1} * [0 i-1 1]';
    [lava_corM{i}, lava_corM_1{i} ,~] = compute_M_M1_synt(lava_flex_info{1}, ortho, 3, tmp(1:3));
    
end


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

% figure;
for ax = 1:size(lava_flex_ax,3)
    
    [x,y,z] = calculate4corners( lava_axM{ax},[0 size2-1], [0 size1-1] );
    
    X_ax_v = [X_ax_v x'];
    Y_ax_v = [Y_ax_v y'];
    Z_ax_v = [Z_ax_v z'];
    
%     points = [x' y' z'];
%     plot3(points(:,1),points(:,2),points(:,3),'b+');hold on
    
    if ax == 1
        N1 = cross([X_ax_v(1,ax) Y_ax_v(1,ax) Z_ax_v(1,ax)]-[X_ax_v(2,ax) Y_ax_v(2,ax) Z_ax_v(2,ax)],[X_ax_v(1,ax) Y_ax_v(1,ax) Z_ax_v(1,ax)]-[X_ax_v(3,ax) Y_ax_v(3,ax) Z_ax_v(3,ax)]); % normal to the axial (ax)
        N1 = N1./norm(N1);
    end
    
end

disp('--------- Calculate M^(-1) and plane eq. for each slice in the third direction -----')

for cor = 1:size1
    
    %     [x,y,z] = calculate4corners_cor( lava_axM{1}, lava_axM{end}, [0 cor-1;511 cor-1]); %  [cor-1 511] , [cor-1 0] [0 cor-1;511 cor-1]
    [x,y,z] = calculate4corners( lava_corM{cor}, [0 size2-1], [0 size3-1]  );
    
    X_cor_v = [X_cor_v x'];
    Y_cor_v = [Y_cor_v y'];
    Z_cor_v = [Z_cor_v z'];
    
%     points = [x' y' z'];
%     plot3(points(:,1),points(:,2),points(:,3),'r*');hold on
    
    if cor == 1
        N3 = cross([X_cor_v(1,cor) Y_cor_v(1,cor) Z_cor_v(1,cor)]-[X_cor_v(2,cor) Y_cor_v(2,cor) Z_cor_v(2,cor)],[X_cor_v(1,cor) Y_cor_v(1,cor) Z_cor_v(1,cor)]-[X_cor_v(3,cor) Y_cor_v(3,cor) Z_cor_v(3,cor)]); % normal to the axial (ax)
        N3 = N3./norm(N3);
    end
    
end

disp('--------- Calculate M^(-1) and plane eq. for each slice in the second direction -----')

for sag = 1:size2
    
    %     [x,y,z] = calculate4corners_cor( lava_axM{1}, lava_axM{end}, [511-sag 511;511-sag 0]); %  [0 511-sag;511 511-sag]
    [x,y,z] = calculate4corners( lava_sagM{sag} , [0 size1-1], [0 size3-1] );
    
    X_sag_v = [X_sag_v x'];
    Y_sag_v = [Y_sag_v y'];
    Z_sag_v = [Z_sag_v z'];
    
%     points = [x' y' z'];
%     plot3(points(:,1),points(:,2),points(:,3),'g*');hold on
    
    if sag == 1
        N2 = cross([X_sag_v(1,sag) Y_sag_v(1,sag) Z_sag_v(1,sag)]-[X_sag_v(2,sag) Y_sag_v(2,sag) Z_sag_v(2,sag)],[X_sag_v(1,sag) Y_sag_v(1,sag) Z_sag_v(1,sag)]-[X_sag_v(3,sag) Y_sag_v(3,sag) Z_sag_v(3,sag)]); % normal to the axial (ax)
        N2 = N2./norm(N2);
    end
    
end

disp('--------- Select the slices for each direction -----')

%% Select a few slices, like in T2W scans
n_slices_ax = 25;
slices_ax = 4:round(size3/n_slices_ax):size3-4;

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

disp('--------- Select the slices axial -----')
%% Apply random deformation to the slices
% figure;
% Axial
for i = 1:length(slices_ax)
    ind = slices_ax(i);
    
    X_ax_def = [X_ax_def X_ax_v(:,ind)];
    Y_ax_def = [Y_ax_def Y_ax_v(:,ind)];
    Z_ax_def = [Z_ax_def Z_ax_v(:,ind)];
    %% image, M, M1
    vol_ax_eval(:,:,i) = lava_flex_ax(:,:,ind);
    axial_M{i}  = lava_axM{ind};
    axial_M1{i} = lava_axM_1{ind};
    
end



% Sagittal
disp('--------- Select the slices sagittal -----')
for i = 1:length(slices_sag)
    ind = slices_sag(i);
    
    X_sag_def = [X_sag_def X_sag_v(:,ind)];
    Y_sag_def = [Y_sag_def Y_sag_v(:,ind)];
    Z_sag_def = [Z_sag_def Z_sag_v(:,ind)];
    %% image, M, M1
    vol_sag_eval(:,:,i) = lava_flex_sag(:,:,ind);
    sag_M{i}  = lava_sagM{ind};
    sag_M1{i} = lava_sagM_1{ind};
    
end


% Coronal
% figure;
disp('--------- Select the slices coronal -----')
for i = 1:length(slices_cor)
    ind = slices_cor(i);
    
    X_cor_def = [X_cor_def X_cor_v(:,ind)];
    Y_cor_def = [Y_cor_def Y_cor_v(:,ind)];
    Z_cor_def = [Z_cor_def Z_cor_v(:,ind)];
    %% image, M, M1
    vol_cor_eval(:,:,i) = lava_flex_cor(:,:,ind);
    
    cor_M{i}  = lava_corM{ind};
    cor_M1{i} = lava_corM_1{ind};
    
end

%% Save the variables in st struct %%
%% 
st.lava_flex_ax  = lava_flex_ax;
st.lava_flex_sag = lava_flex_sag;
st.lava_flex_cor = lava_flex_cor;

%%
st.cor_M  = cor_M;
st.cor_M1 = cor_M1;

st.sag_M  = sag_M;
st.sag_M1 = sag_M1;

st.axial_M  = axial_M;
st.axial_M1 = axial_M1;

%%
st.X_cor = X_cor_def;
st.Y_cor = Y_cor_def;
st.Z_cor = Z_cor_def;

st.X_sag = X_sag_def;
st.Y_sag = Y_sag_def;
st.Z_sag = Z_sag_def;

st.X_ax = X_ax_def;
st.Y_ax = Y_ax_def;
st.Z_ax = Z_ax_def;

%%
st.X_cor_t = X_cor_v;
st.Y_cor_t = Y_cor_v;
st.Z_cor_t = Z_cor_v;

st.X_sag_t = X_sag_v;
st.Y_sag_t = Y_sag_v;
st.Z_sag_t = Z_sag_v;

st.X_ax_t = X_ax_v;
st.Y_ax_t = Y_ax_v;
st.Z_ax_t = Z_ax_v;

st.slices_ax  = slices_ax;
st.slices_sag = slices_sag;
st.slices_cor = slices_cor;

st.vol_ax  = vol_ax_eval;
st.vol_sag = vol_sag_eval;
st.vol_cor = vol_cor_eval;