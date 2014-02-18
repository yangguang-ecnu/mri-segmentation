%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
close all;

target_tri_ax  = TriRep(source_tri_v.Triangulation, xfinal(1:size(source_tri_v.X,1),:));
target_tri_sag = TriRep(source_tri_v.Triangulation, xfinal(size(source_tri_v.X,1)+1:2*size(source_tri_v.X,1),:));
target_tri_cor = TriRep(source_tri_v.Triangulation, xfinal(2*size(source_tri_v.X,1)+1:end,:));

k_ax  = 9;
k_sag = 11;
k_cor = 11;

r_ax = size(vol_ax_eval,1);
c_ax = size(vol_ax_eval,2);

r_sag = size(vol_sag_eval,1);
c_sag = size(vol_sag_eval,2);

r_cor = size(vol_cor_eval,1);
c_cor = size(vol_cor_eval,2);

%% Axial %%

for i = 1:size(vol_ax_eval,1)
    for j = 1:size(vol_ax_eval,2)
        
        p_3d = axial_M{k_ax} * [j i 1]';
        
        current_tr = tsearchn(source_tri_v.X,source_tri_v.Triangulation,[p_3d(1) p_3d(2) p_3d(3)]); % calculate the tetrahedron where p_3d belongs
        
        if ~isnan(current_tr)
            c2b_coord = cartToBary(source_tri_v,current_tr,[p_3d(1) p_3d(2) p_3d(3)]); % get the barycentric coordinates
            b2c_ncoord = baryToCart(target_tri_ax, current_tr, c2b_coord); % get the cartesian coordinates
            
            tmp_v1_ax = axial_M1{k_ax} * [b2c_ncoord 1]'; % 3D point to 2D point in the frame coordinates axial_M_1{sub_3(i)}
            
            %% Make sure the indexes are correct and use bilinear interpolation
            
            fl1ax = floor(tmp_v1_ax(1) + 1);
            fl2ax = floor(tmp_v1_ax(2) + 1);
            cl1ax = ceil(tmp_v1_ax(1)  + 1);
            cl2ax = ceil(tmp_v1_ax(2)  + 1);
            
            min_max_r1ax = min(max(fl2ax,1),r_ax);
            min_max_r2ax = min(max(cl2ax,1),r_ax);
            min_max_c1ax = min(max(fl1ax,1),c_ax);
            min_max_c2ax = min(max(cl1ax,1),c_ax);
            
            neigax = [vol_ax_eval(min_max_r1ax, min_max_c1ax, k_ax) vol_ax_eval(min_max_r1ax, min_max_c2ax, k_ax);...
                      vol_ax_eval(min_max_r2ax, min_max_c1ax, k_ax) vol_ax_eval(min_max_r2ax, min_max_c2ax, k_ax)];
            
            new_axial(i,j,k_ax) = bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax));
        else
            new_axial(i,j,k_ax) = vol_ax_eval(i,j,k_ax);
        end
        
    end
end

% figure;
% subplot(121);imshow(vol_ax_eval(:,:,k),   []);title('Original');
% subplot(122);imshow(new_axial(:,:,k),[]);title('After Deformation');
%% Sagittal %%

for i = 1:size(vol_sag_eval,1)
    for j = 1:size(vol_sag_eval,2)
        
        p_3d = sag_M{k_sag} * [j i 1]';
        
        current_tr = tsearchn(source_tri_v.X,source_tri_v.Triangulation,[p_3d(1) p_3d(2) p_3d(3)]); % calculate the tetrahedron where p_3d belongs
        
        if ~isnan(current_tr)
            c2b_coord = cartToBary(source_tri_v,current_tr,[p_3d(1) p_3d(2) p_3d(3)]); % get the barycentric coordinates
            b2c_ncoord = baryToCart(target_tri_sag, current_tr, c2b_coord); % get the cartesian coordinates
            
            tmp_v1_sag = sag_M1{k_sag} * [b2c_ncoord 1]'; % 3D point to 2D point in the frame coordinates axial_M_1{sub_3(i)}
            
            %% Make sure the indexes are correct and use bilinear interpolation
            
            fl1sag = floor(tmp_v1_sag(1) + 1);
            fl2sag = floor(tmp_v1_sag(2) + 1);
            cl1sag = ceil(tmp_v1_sag(1)  + 1);
            cl2sag = ceil(tmp_v1_sag(2)  + 1);
            
            min_max_r1sag = min(max(fl2sag,1),r_sag);
            min_max_r2sag = min(max(cl2sag,1),r_sag);
            min_max_c1sag = min(max(fl1sag,1),c_sag);
            min_max_c2sag = min(max(cl1sag,1),c_sag);
            
            neigsag = [vol_sag_eval(min_max_r1sag, min_max_c1sag, k_sag) vol_sag_eval(min_max_r1sag, min_max_c2sag, k_sag);...
                       vol_sag_eval(min_max_r2sag, min_max_c1sag, k_sag) vol_sag_eval(min_max_r2sag, min_max_c2sag, k_sag)];
            
            new_sagittal(i,j,k_sag) = bilinear_interpolation(tmp_v1_sag(2), tmp_v1_sag(1), double(neigsag)); % c_sag-j+1
        else
            new_sagittal(i,j,k_sag) = vol_sag_eval(i,j,k_sag); % c_sag-j+1
        end
        
    end
end

% figure;
% subplot(121);imshow(vol_sag_eval(:,:,k),   []);title('Original');
% subplot(122);imshow(new_sagittal(:,:,k),[]);title('After Deformation');

%% Coronal %%

for i = 1:size(vol_cor_eval,1)
    for j = 1:size(vol_cor_eval,2)
        
        p_3d = cor_M{k_cor} * [j i 1]';
        
        current_tr = tsearchn(source_tri_v.X,source_tri_v.Triangulation,[p_3d(1) p_3d(2) p_3d(3)]); % calculate the tetrahedron where p_3d belongs
        
        if ~isnan(current_tr)
            c2b_coord = cartToBary(source_tri_v,current_tr,[p_3d(1) p_3d(2) p_3d(3)]); % get the barycentric coordinates
            b2c_ncoord = baryToCart(target_tri_cor, current_tr, c2b_coord); % get the cartesian coordinates
            
            tmp_v1_cor = cor_M1{k_cor} * [b2c_ncoord 1]'; % 3D point to 2D point in the frame coordinates axial_M_1{sub_3(i)}
            
            %% Make sure the indexes are correct and use bilinear interpolation
            
            fl1cor = floor(tmp_v1_cor(1) + 1);
            fl2cor = floor(tmp_v1_cor(2) + 1);
            cl1cor = ceil(tmp_v1_cor(1)  + 1);
            cl2cor = ceil(tmp_v1_cor(2)  + 1);
            
            min_max_r1cor = min(max(fl2cor,1),r_cor);
            min_max_r2cor = min(max(cl2cor,1),r_cor);
            min_max_c1cor = min(max(fl1cor,1),c_cor);
            min_max_c2cor = min(max(cl1cor,1),c_cor);
            
            neigcor = [vol_cor_eval(min_max_r1cor, min_max_c1cor, k_cor) vol_cor_eval(min_max_r1cor, min_max_c2cor, k_cor);...
                       vol_cor_eval(min_max_r2cor, min_max_c1cor, k_cor) vol_cor_eval(min_max_r2cor, min_max_c2cor, k_cor)];
            
            new_coronal(i,j,k_cor) = bilinear_interpolation(tmp_v1_cor(2), tmp_v1_cor(1), double(neigcor));
        else
            new_coronal(i,j,k_cor) = vol_cor_eval(i,j,k_cor);
        end
        
    end
end

% figure;
% subplot(121);imshow(vol_cor_eval(:,:,k),   []);title('Original');
% subplot(122);imshow(new_coronal(:,:,k),[]);title('After Deformation');

%% Visualize in the RCS
% Original Axial and Sagittal
figure;
aa(:,:,1) = convert2u8(vol_ax_eval(:,:,k_ax));
aa(:,:,2) = convert2u8(vol_ax_eval(:,:,k_ax));
aa(:,:,3) = convert2u8(vol_ax_eval(:,:,k_ax));

[X, Y, Z, triTexture] = compute_RCS_v(axial_M{k_ax},aa,[0 size1-1],[0 size2-1]);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
               'FaceColor','texturemap',...
               'EdgeColor','none');hold on
axis equal

ss(:,:,2) = convert2u8(vol_sag_eval(:,:,k_sag));
ss(:,:,1) = convert2u8(vol_sag_eval(:,:,k_sag));
ss(:,:,3) = convert2u8(vol_sag_eval(:,:,k_sag));

[X, Y, Z, triTexture] = compute_RCS_v(sag_M{k_sag},ss,[0 size3-1],[0 size1-1]);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal

% Original Coronal and Sagittal
figure;
% k_cor = 1;
cc(:,:,2) = convert2u8(vol_cor_eval(:,:,k_cor));
cc(:,:,1) = convert2u8(vol_cor_eval(:,:,k_cor));
cc(:,:,3) = convert2u8(vol_cor_eval(:,:,k_cor));

[X, Y, Z, triTexture] = compute_RCS_v(cor_M{k_cor},cc, [0 size3-1],[0 size2-1]);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal
% k_sag = 24;

ss(:,:,2) = convert2u8(vol_sag_eval(:,:,k_sag));
ss(:,:,1) = convert2u8(vol_sag_eval(:,:,k_sag));
ss(:,:,3) = convert2u8(vol_sag_eval(:,:,k_sag));

[X, Y, Z, triTexture] = compute_RCS_v(sag_M{k_sag},ss,[0 size3-1],[0 size1-1]);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal  

% Original Coronal and Axial
figure;
aa(:,:,1) = convert2u8(vol_ax_eval(:,:,k_ax));
aa(:,:,2) = convert2u8(vol_ax_eval(:,:,k_ax));
aa(:,:,3) = convert2u8(vol_ax_eval(:,:,k_ax));

[X, Y, Z, triTexture] = compute_RCS_v(axial_M{k_ax},aa,[0 size1-1],[0 size2-1]);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
               'FaceColor','texturemap',...
               'EdgeColor','none');hold on
axis equal

cc(:,:,2) = convert2u8(vol_cor_eval(:,:,k_cor));
cc(:,:,1) = convert2u8(vol_cor_eval(:,:,k_cor));
cc(:,:,3) = convert2u8(vol_cor_eval(:,:,k_cor));

[X, Y, Z, triTexture] = compute_RCS_v(cor_M{k_cor},cc, [0 size3-1],[0 size2-1]);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal


% Deformation Axial and Sagittal
figure;
aa(:,:,1) = convert2u8(new_axial(:,:,k_ax));
aa(:,:,2) = convert2u8(new_axial(:,:,k_ax));
aa(:,:,3) = convert2u8(new_axial(:,:,k_ax));

[X, Y, Z, triTexture] = compute_RCS_v(axial_M{k_ax},aa,[0 size1-1],[0 size2-1]);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
               'FaceColor','texturemap',...
               'EdgeColor','none');hold on
axis equal


ss(:,:,2) = convert2u8(new_sagittal(:,:,k_sag));
ss(:,:,1) = convert2u8(new_sagittal(:,:,k_sag));
ss(:,:,3) = convert2u8(new_sagittal(:,:,k_sag));

[X, Y, Z, triTexture] = compute_RCS_v(sag_M{k_sag},ss,[0 size3-1],[0 size1-1]);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal

% Deformation Axial and Coronal

figure
aa(:,:,1) = convert2u8(new_axial(:,:,k_ax));
aa(:,:,2) = convert2u8(new_axial(:,:,k_ax));
aa(:,:,3) = convert2u8(new_axial(:,:,k_ax));

[X, Y, Z, triTexture] = compute_RCS_v(axial_M{k_ax},aa,[0 size1-1],[0 size2-1]);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
               'FaceColor','texturemap',...
               'EdgeColor','none');hold on
axis equal

cc(:,:,2) = convert2u8(new_coronal(:,:,k_cor));
cc(:,:,1) = convert2u8(new_coronal(:,:,k_cor));
cc(:,:,3) = convert2u8(new_coronal(:,:,k_cor));

[X, Y, Z, triTexture] = compute_RCS_v(cor_M{k_cor},cc, [0 size3-1],[0 size2-1]);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal

% Deformation Sagittal and Coronal
figure;

cc(:,:,2) = convert2u8(new_coronal(:,:,k_cor));
cc(:,:,1) = convert2u8(new_coronal(:,:,k_cor));
cc(:,:,3) = convert2u8(new_coronal(:,:,k_cor));

[X, Y, Z, triTexture] = compute_RCS_v(cor_M{k_cor},cc, [0 size3-1],[0 size2-1]);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal


ss(:,:,2) = convert2u8(new_sagittal(:,:,k_sag));
ss(:,:,1) = convert2u8(new_sagittal(:,:,k_sag));
ss(:,:,3) = convert2u8(new_sagittal(:,:,k_sag));

[X, Y, Z, triTexture] = compute_RCS_v(sag_M{k_sag},ss,[0 size3-1],[0 size1-1]);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal
